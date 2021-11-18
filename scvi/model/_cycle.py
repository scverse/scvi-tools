import logging
from typing import Optional

from anndata import AnnData

from scvi.data._anndata import _setup_anndata
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.module import CycleVAE

from .base import BaseModelClass

logger = logging.getLogger(__name__)


class CycleModel(UnsupervisedTrainingMixin, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        **model_kwargs,
    ):
        super().__init__(adata)

        self.module = CycleVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            **model_kwargs,
        )
        self.init_params_ = self._get_init_params(locals())

    @staticmethod
    def setup_anndata(
        adata: AnnData,
        batch_key: Optional[str] = None,
        layer: Optional[str] = None,
        copy: bool = False,
    ) -> Optional[AnnData]:
        return _setup_anndata(
            adata,
            batch_key=batch_key,
            layer=layer,
            copy=copy,
        )
