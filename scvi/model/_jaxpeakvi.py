import logging
from typing import Optional, Sequence

from anndata import AnnData

from scvi.module import JaxPEAKVAE
from scvi.module.base import JaxModuleWrapper
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass, JaxTrainingMixin

logger = logging.getLogger(__name__)



class JaxPEAKVI(JaxTrainingMixin, BaseModelClass):
    
    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        dropout_rate: float = 0.1,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch

        self.module = JaxModuleWrapper(
            JaxPEAKVAE,
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            dropout_rate=dropout_rate,
            **model_kwargs,
        )
        self._model_summary_string = (
            "JaxPEAKVI Model with the following params: \nn_hidden: {}, n_latent: {}, dropout_rate: "
            "{}"
        ).format(
            n_hidden,
            n_latent,
            dropout_rate,
        )
        self.init_params_ = self._get_init_params(locals())
    

    @property
    def device(self):
        return self.module.device
    
    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        **kwargs,
    ):
        pass
    