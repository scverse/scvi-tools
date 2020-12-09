import logging
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from pandas.api.types import CategoricalDtype


from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.data._anndata import _make_obs_column_categorical
from scvi.dataloaders import ConcatDataLoader, ScviDataLoader
from scvi.lightning import SemiSupervisedTask, Trainer, VAETask
from scvi.modules import SCANVAE, VAE

from ._scvi import SCVI
from .base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin
from scvi.lightning._sampling import SubSampleLabels

logger = logging.getLogger(__name__)


class SCANVI(RNASeqMixin, VAEMixin, ArchesMixin, BaseModelClass):
    """
    Single-cell annotation using variational inference [Xu20]_.

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
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.modules.VAE` and :class:`~scvi.modules.SCANVAE`

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
        use_gpu: bool = True,
        vae_model_kwargs: dict = {},
        scanvae_model_kwargs: dict = {},
    ):
        super(SCANVI, self).__init__(adata, use_gpu=use_gpu)
        vae_model_kwargs = dict(vae_model_kwargs)
        scanvae_model_kwargs = dict(scanvae_model_kwargs)
        self.unlabeled_category_ = unlabeled_category
        has_unlabeled = self._set_indices_and_labels()

        if pretrained_model is not None:
            if pretrained_model.is_trained is False:
                raise ValueError("pretrained model has not been trained")
            self._base_model = pretrained_model.model
            self._is_trained_base = True
        else:
            self._base_model = VAE(
                n_input=self.summary_stats["n_vars"],
                n_batch=self.summary_stats["n_batch"],
                n_hidden=n_hidden,
                n_latent=n_latent,
                n_layers=n_layers,
                dropout_rate=dropout_rate,
                dispersion=dispersion,
                gene_likelihood=gene_likelihood,
                **vae_model_kwargs,
            )
            self._is_trained_base = False

        # ignores unlabelled catgegory
        n_labels = (
            self.summary_stats["n_labels"] - 1
            if has_unlabeled
            else self.summary_stats["n_labels"]
        )
        self.model = SCANVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            **scanvae_model_kwargs,
        )

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            "ScanVI Model with the following params: \nunlabeled_category: {}, n_hidden: {}, n_latent: {}"
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

    def _set_indices_and_labels(self):
        """
        Set indices and make unlabeled cat as the last cat.

        Returns
        -------
        True is categories reordered else False
        """
        # get indices for labeled and unlabeled cells
        key = self.scvi_setup_dict_["data_registry"][_CONSTANTS.LABELS_KEY]["attr_key"]
        mapping = self.scvi_setup_dict_["categorical_mappings"][key]["mapping"]
        original_key = self.scvi_setup_dict_["categorical_mappings"][key][
            "original_key"
        ]
        labels = np.asarray(self.adata.obs[original_key]).ravel()

        if self.unlabeled_category_ in labels:
            unlabeled_idx = np.where(mapping == self.unlabeled_category_)
            unlabeled_idx = unlabeled_idx[0][0]
            # move unlabelled category to be the last position
            mapping[unlabeled_idx], mapping[-1] = mapping[-1], mapping[unlabeled_idx]
            cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
            # rerun setup for the batch column
            _make_obs_column_categorical(
                self.adata,
                original_key,
                "_scvi_labels",
                categorical_dtype=cat_dtype,
            )
            remapped = True
        else:
            remapped = False

        self.scvi_setup_dict_ = self.adata.uns["_scvi"]
        self._label_mapping = mapping
        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(
            labels == self.unlabeled_category_
        ).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category_).ravel()
        self._code_to_label = {i: l for i, l in enumerate(self._label_mapping)}
        return remapped

    @property
    def _task_class(self):
        return SemiSupervisedTask

    @property
    def _data_loader_cls(self):
        return ScviDataLoader

    @property
    def history(self):
        """Returns computed metrics during training."""
        return {
            "unsupervised_trainer_history": self.unsupervised_history_,
            "semisupervised_trainer_history": self.semisupervised_history_,
        }

    def train(
        self,
        n_epochs_unsupervised: Optional[int] = None,
        n_epochs_semisupervised: Optional[int] = None,
        n_samples_per_label: Optional[float] = None,
        validation_size: Optional[float] = None,
        batch_size: int = 512,
        use_gpu: Optional[bool] = None,
        num_workers: int = 4,
        train_base_model: bool = True,
        lr: float = 1e-3,
        n_epochs_kl_warmup: int = 400,
        n_iter_kl_warmup: Optional[int] = None,
        check_val_every_n_epoch: Optional[int] = None,
        unsupervised_trainer_kwargs: dict = {},
        semisupervised_trainer_kwargs: dict = {},
        unsupervised_task_kwargs: dict = {},
        semisupervised_task_kwargs: dict = {},
    ):
        """
        Train the model.

        Parameters
        ----------
        n_epochs_unsupervised
            Number of passes through the dataset for unsupervised pre-training.
        n_epochs_semisupervised
            Number of passes through the dataset for semisupervised training.
        train_base_model
            Pretrain an SCVI base model first before semisupervised training.
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
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        unsupervised_trainer_kwargs
            Other keyword args for :class:`~scvi.trainers.UnsupervisedTrainer`.
        semisupervised_trainer_kwargs
            Other keyword args for :class:`~scvi.trainers.SemiSupervisedTrainer`.
        semisupervised_train_kwargs
            Keyword args for the train method of :class:`~scvi.trainers.SemiSupervisedTrainer`.
        """
        unsupervised_trainer_kwargs = dict(unsupervised_trainer_kwargs)
        semisupervised_trainer_kwargs = dict(semisupervised_trainer_kwargs)
        unsupervised_task_kwargs = dict(unsupervised_task_kwargs)
        semisupervised_task_kwargs = dict(semisupervised_task_kwargs)

        if n_epochs_unsupervised is None:
            n_epochs_unsupervised = np.min(
                [round((20000 / self.adata.shape[0]) * 400), 400]
            )
        if n_epochs_semisupervised is None:
            n_epochs_semisupervised = int(
                np.min([10, np.max([2, round(n_epochs_unsupervised / 3.0)])])
            )
        logger.info(
            "Training Unsupervised Trainer for {} epochs.".format(n_epochs_unsupervised)
        )
        logger.info(
            "Training SemiSupervised Trainer for {} epochs.".format(
                n_epochs_semisupervised
            )
        )
        if use_gpu is not None:
            use_gpu = use_gpu
        else:
            use_gpu = self.use_gpu

        if use_gpu:
            gpus = 1
            pin_memory = True
        else:
            gpus = None
            pin_memory = False
        self.trainer = Trainer(
            max_epochs=n_epochs_unsupervised,
            gpus=gpus,
            **unsupervised_trainer_kwargs,
        )

        # self.test_indices_ = test_dl.indices
        # self.validation_indices_ = val_dl.indices
        full_idx = np.arange(self.adata.n_obs)
        train_dl = ScviDataLoader(
            self.adata,
            shuffle=True,
            indices=full_idx,
            batch_size=batch_size,
            pin_memory=pin_memory,
            num_workers=num_workers,
        )
        self.train_indices_ = train_dl.indices
        self._unsupervised_task = VAETask(
            self._base_model, len(self.train_indices), **unsupervised_task_kwargs
        )
        self.trainer.fit(self._unsupervised_task, train_dl)
        self._base_model.eval()

        self.model.load_state_dict(self._base_model.state_dict(), strict=False)

        self._semisupervised_task = SemiSupervisedTask(self.model)

        # if we have labelled cells, we want to pass them through the classifier
        # else we only pass the full dataset
        if len(self._labeled_indices) != 0:
            semisupervised_train_dl = SemiSupervisedDataLoader(
                self.adata,
                full_idx,
                self.original_label_key,
                self.unlabeled_category_,
                n_samples_per_label,
                shuffle=True,
                batch_size=batch_size,
                pin_memory=pin_memory,
                num_workers=num_workers,
            )
            sampler_callback = [SubSampleLabels()]
        else:
            semisupervised_train_dl = ScviDataLoader(
                self.adata,
                shuffle=True,
                indices=full_idx,
                batch_size=batch_size,
                pin_memory=pin_memory,
                num_workers=num_workers,
            )
            sampler_callback = None
        self._semisupervised_trainer = Trainer(
            max_epochs=n_epochs_semisupervised,
            gpus=gpus,
            callbacks=sampler_callback,
        )

        self._semisupervised_trainer.fit(
            self._semisupervised_task, semisupervised_train_dl
        )
        self.model.eval()
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
        y_pred = []
        for _, tensors in enumerate(scdl):
            x = tensors[_CONSTANTS.X_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            pred = self.model.classify(x, batch)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

        y_pred = np.array(torch.cat(y_pred))
        if not soft:
            predictions = []
            for p in y_pred:
                predictions.append(self._code_to_label[p])

            return np.array(predictions)
        else:
            n_labels = len(pred[0])
            pred = pd.DataFrame(
                pred,
                columns=self._label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
            return y_pred
