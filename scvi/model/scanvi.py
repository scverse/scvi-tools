import logging
import numpy as np
import pandas as pd

from typing import Union, Optional, Sequence
from anndata import AnnData

from scvi import _CONSTANTS
from scvi.model import SCVI
from scvi.core.models import BaseModelClass, RNASeqMixin, VAEMixin
from scvi.core.modules import VAE, SCANVAE
from scvi.core.trainers import UnsupervisedTrainer, SemiSupervisedTrainer
from scvi.core.data_loaders import AnnotationDataLoader
from scvi._compat import Literal

logger = logging.getLogger(__name__)


class SCANVI(RNASeqMixin, VAEMixin, BaseModelClass):
    """
    Single-cell annotation using variational inference [Xu19]_.

    Inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    unlabeled_category
        Value used for unlabeled cells in `labels_key` used to setup AnnData with scvi.
    pretrained_model
        Instance of SCVI model that has already been trained.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    use_cuda
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.core.modules.VAE` and :class:`~scvi.core.modules.SCANVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.data.setup_anndata(adata, batch_key="batch", labels_key="labels")
    >>> vae = scvi.model.SCANVI(adata, "Unknown")
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obs["pred_label"] = vae.predict()
    """

    def __init__(
        self,
        adata: AnnData,
        unlabeled_category: Union[str, int, float],
        pretrained_model: Optional[SCVI] = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        use_cuda: bool = True,
        **model_kwargs,
    ):
        super(SCANVI, self).__init__(adata, use_cuda=use_cuda)
        self.unlabeled_category = unlabeled_category

        if pretrained_model is not None:
            if pretrained_model.is_trained is False:
                raise ValueError("pretrained model has not been trained")
            self._base_model = pretrained_model.model
            self._is_trained_base = True
        else:
            self._base_model = VAE(
                n_input=self.summary_stats["n_genes"],
                n_batch=self.summary_stats["n_batch"],
                n_hidden=n_hidden,
                n_latent=n_latent,
                n_layers=n_layers,
                dropout_rate=dropout_rate,
                dispersion=dispersion,
                gene_likelihood=gene_likelihood,
                **model_kwargs,
            )
            self._is_trained_base = False
        self.model = SCANVAE(
            n_input=self.summary_stats["n_genes"],
            n_batch=self.summary_stats["n_batch"],
            n_labels=self.summary_stats["n_labels"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            **model_kwargs,
        )

        # get indices for labeled and unlabeled cells
        key = self.scvi_setup_dict_["data_registry"][_CONSTANTS.LABELS_KEY]["attr_key"]
        self._label_mapping = self.scvi_setup_dict_["categorical_mappings"][key][
            "mapping"
        ]
        original_key = self.scvi_setup_dict_["categorical_mappings"][key][
            "original_key"
        ]
        labels = np.asarray(self.adata.obs[original_key]).ravel()
        self._code_to_label = {i: l for i, l in enumerate(self._label_mapping)}
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category).ravel()

        self._model_summary_string = (
            "ScanVI Model with params: \nunlabeled_category: {}, n_hidden: {}, n_latent: {}"
            ", n_layers: {}, dropout_rate: {}, dispersion: {}, gene_likelihood: {}"
        ).format(
            unlabeled_category,
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            gene_likelihood,
        )
        self.init_params_ = self._get_init_params(locals())

    @property
    def _trainer_class(self):
        return SemiSupervisedTrainer

    @property
    def _scvi_dl_class(self):
        return AnnotationDataLoader

    def train(
        self,
        n_epochs_unsupervised=None,
        n_epochs_semisupervised=None,
        train_size=0.9,
        test_size=None,
        lr=1e-3,
        n_epochs_kl_warmup=400,
        n_iter_kl_warmup=None,
        frequency=None,
        unsupervised_trainer_kwargs={},
        semisupervised_trainer_kwargs={},
        unsupervised_train_kwargs={},
        semisupervised_train_kwargs={},
    ):
        """
        Train the model.

        Parameters
        ----------
        n_epochs_unsupervised
            Number of passes through the dataset for unsupervised pre-training.
        n_epochs_semisupervised
            Number of passes through the dataset for semisupervised training.
        train_size
            Size of training set in the range [0.0, 1.0].
        test_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + test_size < 1`, the remaining cells belong to a validation set.
        lr
            Learning rate for optimization.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        n_iter_kl_warmup
            Number of minibatches for scaling term on KL divergence to go from 0 to 1.
            To use, set to not `None` and set `n_epochs_kl_warmup` to `None`.
        frequency
            Frequency with which metrics are computed on the data for train/test/val sets.
        unsupervised_trainer_kwargs
            Other keyword args for :class:`~scvi.core.trainers.UnsupervisedTrainer`.
        semisupervised_trainer_kwargs
            Other keyword args for :class:`~scvi.core.trainers.SemiSupervisedTrainer`.
        semisupervised_train_kwargs
            Keyword args for the train method of :class:`~scvi.core.trainers.SemiSupervisedTrainer`.
        """
        unsupervised_trainer_kwargs = dict(unsupervised_trainer_kwargs)
        semisupervised_trainer_kwargs = dict(semisupervised_trainer_kwargs)
        unsupervised_train_kwargs = dict(unsupervised_train_kwargs)
        semisupervised_train_kwargs = dict(semisupervised_train_kwargs)

        if n_epochs_unsupervised is None:
            n_epochs_unsupervised = np.min(
                [round((20000 / self.adata.shape[0]) * 400), 400]
            )
        if n_epochs_semisupervised is None:
            n_epochs_semisupervised = int(
                np.min([10, np.max([2, round(n_epochs_unsupervised / 3.0)])])
            )

        if self._is_trained_base is not True:
            self._unsupervised_trainer = UnsupervisedTrainer(
                self._base_model,
                self.adata,
                train_size=train_size,
                test_size=test_size,
                n_iter_kl_warmup=n_iter_kl_warmup,
                n_epochs_kl_warmup=n_epochs_kl_warmup,
                frequency=frequency,
                use_cuda=self.use_cuda,
                **unsupervised_trainer_kwargs,
            )
            self._unsupervised_trainer.train(
                n_epochs=n_epochs_unsupervised, lr=lr, **unsupervised_train_kwargs
            )
            self._is_trained_base = True

        self.model.load_state_dict(self._base_model.state_dict(), strict=False)

        self.trainer = SemiSupervisedTrainer(
            self.model,
            self.adata,
            use_cuda=self.use_cuda,
            **semisupervised_trainer_kwargs,
        )

        self.trainer.unlabelled_set = self.trainer.create_scvi_dl(
            indices=self._unlabeled_indices
        )
        self.trainer.labelled_set = self.trainer.create_scvi_dl(
            indices=self._labeled_indices
        )
        self.trainer.train(
            n_epochs=n_epochs_semisupervised, **semisupervised_train_kwargs
        )

        self.is_trained_ = True

    def predict(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        soft: bool = False,
        batch_size: int = 128,
    ) -> Union[np.ndarray, pd.DataFrame]:
        """
        Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        soft
            Return probabilities for each class label.
        batch_size
            Minibatch size to use.
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        _, pred = scdl.sequential().compute_predictions(soft=soft)

        if not soft:
            predictions = []
            for p in pred:
                predictions.append(self._code_to_label[p])

            return np.array(predictions)
        else:
            pred = pd.DataFrame(
                pred,
                columns=self._label_mapping,
                index=adata.obs_names[indices],
            )
            return pred
