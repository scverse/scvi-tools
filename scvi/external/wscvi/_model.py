import logging

from anndata import AnnData

from scvi._compat import Literal
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.external.wscvi._module import WVAE

from scvi.model.base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin

logger = logging.getLogger(__name__)


class WSCVI(
    RNASeqMixin, VAEMixin, ArchesMixin, UnsupervisedTrainingMixin, BaseModelClass
):
    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super(WSCVI, self).__init__(adata)

        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )
        self.module = WVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
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
            gene_likelihood,
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())
