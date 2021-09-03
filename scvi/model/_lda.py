import logging
from typing import Optional

import numpy as np
from anndata import AnnData

from scvi.module import LDAPyroModule

from .base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class LDA(PyroSviTrainMixin, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        n_components: int = 20,
        n_hidden: int = 128,
        cell_component_prior: Optional[np.ndarray] = None,
        component_gene_prior: Optional[np.ndarray] = None,
    ):
        super().__init__(adata)

        self.module = LDAPyroModule(
            n_input=self.summary_stats["n_vars"],
            n_components=n_components,
            n_hidden=n_hidden,
            cell_component_prior=cell_component_prior,
            component_gene_prior=component_gene_prior,
        )
