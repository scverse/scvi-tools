import logging

from anndata import AnnData

from scvi.module import LDAModule

from .base import BaseModelClass

logger = logging.getLogger(__name__)


class LDA(BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        n_topics: int = 20,
        n_hidden: int = 128,
    ):
        super().__init__(adata)

        self.module = LDAModule(
            n_input=self.summary_stats["n_vars"],
            n_topics=n_topics,
            n_hidden=n_hidden,
        )
