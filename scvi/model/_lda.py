import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData

from scvi._constants import _CONSTANTS
from scvi.dataloaders._anntorchdataset import AnnTorchDataset
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
        # in case any other model was created before that shares the same parameter names.
        pyro.clear_param_store()

        super().__init__(adata)

        self.module = LDAPyroModule(
            n_input=self.summary_stats["n_vars"],
            n_components=n_components,
            n_hidden=n_hidden,
            cell_component_prior=cell_component_prior,
            component_gene_prior=component_gene_prior,
        )

    def _check_var_equality(self, adata: AnnData):
        source_var_names = self.adata.var_names.astype(str)
        user_var_names = adata.var_names.astype(str)
        if not np.array_equal(user_var_names, source_var_names):
            raise ValueError(
                "`adata` passed into `transform` does not have matching var_names "
                "with the source adata the model was trained with."
            )

    def get_components(self) -> pd.DataFrame:
        if self.is_trained_ is False:
            warnings.warn(
                "Trying to query inferred values from an untrained model. Please train the model first."
            )

        return pd.DataFrame(
            data=self.module.components.numpy(), index=self.adata.obs_names
        )

    def transform(self, adata: Optional[AnnData] = None) -> pd.DataFrame:
        if adata is not None:
            self._check_var_equality(adata)
        user_adata = adata or self.adata
        # TODO(jhong): remove jankiness
        adata_dataset = AnnTorchDataset(user_adata)

        transformed_X = self.module.transform(
            torch.from_numpy(adata_dataset.get_data(_CONSTANTS.X_KEY))
        ).numpy()
        return pd.DataFrame(data=transformed_X, index=user_adata.obs_names)

    def perplexity(self, adata: Optional[AnnData] = None) -> float:
        if adata is not None:
            self._check_var_equality(adata)
        user_adata = adata or self.adata
        adata_dataset = AnnTorchDataset(user_adata)

        return self.module.perplexity(
            torch.from_numpy(adata_dataset.get_data(_CONSTANTS.X_KEY))
        )
