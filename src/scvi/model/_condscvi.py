from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
import torch

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalObsField,
)
from scvi.model.base import (
    BaseMinifiedModeModelClass,
    RNASeqMixin,
    SemisupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module import VAEC
from scvi.utils import setup_anndata_dsp

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, SemisupervisedTrainingMixin, BaseMinifiedModeModelClass):
    """Conditional version of single-cell Variational Inference.

    Used for multi-resolution deconvolution of spatial transcriptomics data :cite:p:`Lopez22`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.CondSCVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    weight_obs
        Whether to reweight observations by their inverse proportion (useful for lowly abundant
        cell types)
    dropout_rate
        Dropout rate for neural networks.
    **module_kwargs
        Keyword args for :class:`~scvi.module.VAEC`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.CondSCVI.setup_anndata(adata, "labels")
    >>> vae = scvi.model.CondSCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_CondSCVI"] = vae.get_latent_representation()

    Notes
    -----
    See further usage examples in the following tutorial:

    1. :doc:`/tutorials/notebooks/spatial/DestVI_tutorial`
    """

    _module_cls = VAEC
    _LATENT_QZM_KEY = "condscvi_latent_qzm"
    _LATENT_QZV_KEY = "condscvi_latent_qzv"

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        weight_obs: bool = False,
        dropout_rate: float = 0.05,
        **module_kwargs,
    ):
        super().__init__(adata)

        self.n_labels = self.summary_stats.n_labels
        self.n_vars = self.summary_stats.n_vars
        self._set_indices_and_labels()
        if weight_obs:
            ct_counts = np.unique(
                self.get_from_registry(adata, REGISTRY_KEYS.LABELS_KEY),
                return_counts=True,
            )[1]
            ct_prop = ct_counts / np.sum(ct_counts)
            ct_prop[ct_prop < 0.05] = 0.05
            ct_prop = ct_prop / np.sum(ct_prop)
            ct_weight = 1.0 / ct_prop
            module_kwargs.update({"ct_weight": ct_weight})

        self.module = self._module_cls(
            n_input=self.n_vars,
            n_batch=getattr(self.summary_stats, "n_batch", 0),
            n_labels=self.n_labels,
            n_fine_labels=getattr(self.summary_stats, "n_fine_labels", None),
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **module_kwargs,
        )
        self._model_summary_string = (
            f"Conditional SCVI Model with the following params: \nn_hidden: {n_hidden}, "
            f"n_latent: {n_latent}, n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"weight_obs: {weight_obs}"
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def get_vamp_prior(self, adata: AnnData | None = None, p: int = 10) -> np.ndarray:
        r"""Return an empirical prior over the cell-type specific latent space (vamp prior).

        May be used for deconvolution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        p
            number of clusters in kmeans clustering for cell-type sub-clustering for empirical
            prior

        Returns
        -------
        mean_vprior: np.ndarray
            (n_labels, p, D) array
        var_vprior
            (n_labels, p, D) array
        weights_vprior
            (n_labels, p) array
        """
        from sklearn.cluster import KMeans

        if self.is_trained_ is False:
            warnings.warn(
                "Trying to query inferred values from an untrained model. "
                "Please train the model first.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        adata = self._validate_anndata(adata)

        if self.module.prior == "mog":
            results = {
                "mean_vprior": self.module.prior_means,
                "var_vprior": torch.exp(self.module.prior_log_std) ** 2,
                "weights_vprior": torch.nn.functional.softmax(self.module.prior_logits, dim=-1),
            }
        else:
            # Extracting latent representation of adata including variances.
            mean_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
            var_vprior = np.ones((self.summary_stats.n_labels, p, self.module.n_latent))
            mp_vprior = np.zeros((self.summary_stats.n_labels, p))

            labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
            key = labels_state_registry.original_key
            mapping = labels_state_registry.categorical_mapping

            mean_cat, var_cat = self.get_latent_representation(
                adata, return_dist=True, batch_size=p, give_mean=False
            )

            for ct in range(self.summary_stats["n_labels"]):
                local_indices = np.where(adata.obs[key] == mapping[ct])[0]
                n_local_indices = len(local_indices)
                if "overclustering_vamp" not in adata.obs.columns:
                    if p < n_local_indices and p > 0:
                        overclustering_vamp = KMeans(n_clusters=p, n_init=30).fit_predict(
                            mean_cat[local_indices]
                        )
                    else:
                        # Every cell is its own cluster
                        overclustering_vamp = np.arange(n_local_indices)
                else:
                    overclustering_vamp = adata[local_indices, :].obs["overclustering_vamp"]

                keys, counts = np.unique(overclustering_vamp, return_counts=True)

                n_labels_overclustering = len(keys)
                if n_labels_overclustering > p:
                    error_mess = """
                        Given cell type specific clustering contains more clusters than
                        vamp_prior_p. Increase value of vamp_prior_p to largest number of cell type
                         specific clusters."""

                    raise ValueError(error_mess)

                var_cluster = np.ones(
                    [
                        n_labels_overclustering,
                        self.module.n_latent,
                    ]
                )
                mean_cluster = np.zeros_like(var_cluster)

                for index, cluster in enumerate(keys):
                    indices_curr = local_indices[np.where(overclustering_vamp == cluster)[0]]
                    var_cluster[index, :] = np.mean(var_cat[indices_curr], axis=0) + np.var(
                        mean_cat[indices_curr], axis=0
                    )
                    mean_cluster[index, :] = np.mean(mean_cat[indices_curr], axis=0)

                slicing = slice(n_labels_overclustering)
                mean_vprior[ct, slicing, :] = mean_cluster
                var_vprior[ct, slicing, :] = var_cluster
                mp_vprior[ct, slicing] = counts / sum(counts)
            results = {
                "mean_vprior": mean_vprior,
                "var_vprior": var_vprior,
                "weights_vprior": mp_vprior,
            }

        return results

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        fine_labels_key: str | None = None,
        layer: str | None = None,
        unlabeled_category: str = "unlabeled",
        size_factor_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        %(param_labels_key)s
        fine_labels_key
            Key in `adata.obs` where fine-grained labels are stored.
        %(param_layer)s
        %(param_unlabeled_category)s
        %(param_size_factor_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
        ]
        if fine_labels_key is not None:
            anndata_fields.append(
                LabelsWithUnlabeledObsField("fine_labels", fine_labels_key, unlabeled_category),
            )

        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(cls, adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
