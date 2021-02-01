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
from scvi.dataloaders import AnnDataLoader, SemiSupervisedDataLoader
from scvi.modules import SCANVAE
from .base import (
    ArchesMixin,
    BaseModelClass,
    RNASeqMixin,
    VAEMixin,
    SemiSupervisedTrainingMixin,
)
from .base._utils import _load_torch_weights

logger = logging.getLogger(__name__)


class SCANVI(
    RNASeqMixin, VAEMixin, ArchesMixin, SemiSupervisedTrainingMixin, BaseModelClass
):
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
    def history(self):
        """Returns computed metrics during training."""
        return self._trainer.logger.history

    def load_pretrained_scvi(self, dir_path: str):
        map_location = torch.device("cpu") if self.use_gpu is False else None
        scvi_state_dict = _load_torch_weights(dir_path, map_location)
        self.module.load_state_dict(scvi_state_dict, strict=False)

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
