import logging
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from pandas.api.types import CategoricalDtype
from sklearn.model_selection._split import _validate_shuffle_split

from scvi import _CONSTANTS, settings
from scvi._compat import Literal
from scvi.data._anndata import _make_obs_column_categorical
from scvi.dataloaders import AnnDataLoader, SemiSupervisedDataLoader
from scvi.lightning import SemiSupervisedTrainingPlan, Trainer
from scvi.lightning._callbacks import SubSampleLabels
from scvi.modules import SCANVAE

from .base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin

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
        Keyword args for :class:`~scvi.modules.SCANVAE`

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
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(SCANVI, self).__init__(adata, use_gpu=use_gpu)
        scanvae_model_kwargs = dict(model_kwargs)

        self.unlabeled_category_ = unlabeled_category
        has_unlabeled = self._set_indices_and_labels()

        if len(self._labeled_indices) != 0:
            self._dl_cls = SemiSupervisedDataLoader
        else:
            self._dl_cls = AnnDataLoader

        # ignores unlabeled catgegory
        n_labels = (
            self.summary_stats["n_labels"] - 1
            if has_unlabeled
            else self.summary_stats["n_labels"]
        )
        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )
        self.module = SCANVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_labels=n_labels,
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
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
            # move unlabeled category to be the last position
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
        self.original_label_key = original_key

        return remapped

    @property
    def _plan_class(self):
        return SemiSupervisedTrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    @property
    def history(self):
        """Returns computed metrics during training."""
        return self._trainer.logger.history

    def train(
        self,
        max_epochs: Optional[int] = None,
        n_samples_per_label: Optional[float] = None,
        check_val_every_n_epoch: Optional[int] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        use_gpu: Optional[bool] = None,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        use_gpu
            If `True`, use the GPU if available. Will override the use_gpu option when initializing model
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.SemiSupervisedTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        logger.info("Training for {} epochs.".format(max_epochs))

        use_gpu = use_gpu if use_gpu is not None else self.use_gpu
        gpus = 1 if use_gpu else None
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and use_gpu) else False
        )
        train_dl, val_dl, test_dl = self._train_test_val_split(
            self.adata,
            unlabeled_category=self.unlabeled_category_,
            train_size=train_size,
            validation_size=validation_size,
            n_samples_per_label=n_samples_per_label,
            pin_memory=pin_memory,
            batch_size=batch_size,
        )

        self.train_indices_ = train_dl.indices
        self.validation_indices_ = val_dl.indices
        self.test_indices_ = test_dl.indices

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs
        self._task = SemiSupervisedTrainingPlan(self.module, **plan_kwargs)

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = (
            [SubSampleLabels()] if len(self._labeled_indices) != 0 else []
        )
        self._trainer = Trainer(
            max_epochs=max_epochs,
            gpus=gpus,
            callbacks=sampler_callback,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **kwargs,
        )
        if len(self.validation_indices_) != 0:
            self._trainer.fit(self._task, train_dl, val_dl)
        else:
            self._trainer.fit(self._task, train_dl)
        self.module.eval()
        self.is_trained_ = True

    def predict(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        soft: bool = False,
        batch_size: Optional[int] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        """
        Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
        indices
            Return probabilities for each class label.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_scvi_dl(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )
        y_pred = []
        for _, tensors in enumerate(scdl):
            x = tensors[_CONSTANTS.X_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            pred = self.module.classify(x, batch)
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
                y_pred,
                columns=self._label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
            return y_pred

    def _train_test_val_split(
        self,
        adata: AnnData,
        unlabeled_category,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        n_samples_per_label: Optional[int] = None,
        **kwargs,
    ):
        """
        Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

        If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.
        The ratio between labeled and unlabeled data in adata will be preserved
        in the train/test/val sets.

        Parameters
        ----------
        adata
            AnnData to split into train/test/val sets
        unlabeled_category
            Category to treat as unlabeled
        train_size
            float, or None (default is 0.9)
        validation_size
            float, or None (default is None)
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch
        **kwargs
            Keyword args for `_make_scvi_dl()`
        """
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        n_labeled_idx = len(self._labeled_indices)
        n_unlabeled_idx = len(self._unlabeled_indices)

        def get_train_val_split(n_samples, test_size, train_size):
            try:
                n_train, n_val = _validate_shuffle_split(
                    n_samples, test_size, train_size
                )
            except ValueError:
                if train_size != 1.0 and n_samples != 1:
                    raise ValueError(
                        "Choice of train_size={} and validation_size={} not understood".format(
                            train_size, test_size
                        )
                    )
                n_train, n_val = n_samples, 0
            return n_train, n_val

        if n_labeled_idx != 0:
            n_labeled_train, n_labeled_val = get_train_val_split(
                n_labeled_idx, validation_size, train_size
            )
            labeled_permutation = np.random.choice(
                self._labeled_indices, len(self._labeled_indices), replace=False
            )
            labeled_idx_val = labeled_permutation[:n_labeled_val]
            labeled_idx_train = labeled_permutation[
                n_labeled_val : (n_labeled_val + n_labeled_train)
            ]
            labeled_idx_test = labeled_permutation[(n_labeled_val + n_labeled_train) :]
        else:
            labeled_idx_test = []
            labeled_idx_train = []
            labeled_idx_val = []

        if n_unlabeled_idx != 0:
            n_unlabeled_train, n_unlabeled_val = get_train_val_split(
                n_unlabeled_idx, validation_size, train_size
            )
            unlabeled_permutation = np.random.choice(
                self._unlabeled_indices, len(self._unlabeled_indices)
            )
            unlabeled_idx_val = unlabeled_permutation[:n_unlabeled_val]
            unlabeled_idx_train = unlabeled_permutation[
                n_unlabeled_val : (n_unlabeled_val + n_unlabeled_train)
            ]
            unlabeled_idx_test = unlabeled_permutation[
                (n_unlabeled_val + n_unlabeled_train) :
            ]
        else:
            unlabeled_idx_train = []
            unlabeled_idx_val = []
            unlabeled_idx_test = []

        indices_train = np.concatenate((labeled_idx_train, unlabeled_idx_train))
        indices_val = np.concatenate((labeled_idx_val, unlabeled_idx_val))
        indices_test = np.concatenate((labeled_idx_test, unlabeled_idx_test))

        indices_train = indices_train.astype(int)
        indices_val = indices_val.astype(int)
        indices_test = indices_test.astype(int)

        if len(self._labeled_indices) != 0:
            dataloader_class = SemiSupervisedDataLoader
            dl_kwargs = {
                "unlabeled_category": unlabeled_category,
                "n_samples_per_label": n_samples_per_label,
            }
        else:
            dataloader_class = AnnDataLoader
            dl_kwargs = {}
        dl_kwargs.update(kwargs)

        scanvi_train_dl = self._make_scvi_dl(
            adata,
            indices=indices_train,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )
        scanvi_val_dl = self._make_scvi_dl(
            adata,
            indices=indices_val,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )
        scanvi_test_dl = self._make_scvi_dl(
            adata,
            indices=indices_test,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )

        return scanvi_train_dl, scanvi_val_dl, scanvi_test_dl
