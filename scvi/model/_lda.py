import logging

from anndata import AnnData

from scvi.module import LDAPyroModule

from .base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class LDA(PyroSviTrainMixin, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        n_topics: int = 20,
        n_hidden: int = 128,
    ):
        super().__init__(adata)

        self.module = LDAPyroModule(
            n_input=self.summary_stats["n_vars"],
            n_topics=n_topics,
            n_hidden=n_hidden,
        )
