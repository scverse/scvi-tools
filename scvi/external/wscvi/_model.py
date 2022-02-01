import logging
from typing import List, Optional

from anndata import AnnData

from scvi._compat import Literal
from scvi.data._anndata import _setup_anndata
from scvi.external.wscvi._module import WVAE
from scvi.model._utils import _init_library_size
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    DEMixin,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)


class WSCVI(
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    DEMixin,
    BaseModelClass,
):
    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.0,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super(WSCVI, self).__init__(adata)
        self.n_genes = adata.shape[-1]
        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )
        n_batch = self.summary_stats["n_batch"]
        library_log_means, library_log_vars = _init_library_size(adata, n_batch)

        self.module = WVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=n_batch,
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            latent_distribution=latent_distribution,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            **model_kwargs,
        )
        self._model_summary_string = (
            "SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            "nb",
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    @staticmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        adata: AnnData,
        batch_key: Optional[str] = None,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        copy: bool = False,
    ) -> Optional[AnnData]:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_layer)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        %(param_copy)s

        Returns
        -------
        %(returns)s
        """
        return _setup_anndata(
            adata,
            batch_key=batch_key,
            labels_key=labels_key,
            layer=layer,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys,
            copy=copy,
        )
