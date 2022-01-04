import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ._obsm_field import ObsmField

logger = logging.getLogger(__name__)


class ProteinObsmField(ObsmField):
    """
    An AnnDataField for an protein data stored in an .obsm field of an AnnData object.

    For usage with the TotalVI model. Computes an additional mask which indicates
    where batches are missing protein data.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obsm_key
        Key to access the field in the AnnData .obsm mapping.
    batch_key
        Key corresponding to the .obs field where batch indices are stored.
        Used for computing a batch mask over the data for missing protein data.
    colnames_uns_key
        Key to access column names corresponding to each column of the .obsm field in
        the AnnData .uns mapping.
    is_count_data
        If True, checks if the data are counts during validation.
    """

    PROTEIN_BATCH_MASK = "protein_batch_mask"

    def __init__(
        self,
        registry_key: str,
        obsm_key: str,
        batch_key: str,
        colnames_uns_key: Optional[str] = None,
        is_count_data: bool = False,
    ) -> None:
        super().__init__(
            registry_key,
            obsm_key,
            colnames_uns_key=colnames_uns_key,
            is_count_data=is_count_data,
        )
        self.batch_key = batch_key

    def _get_batch_mask_protein_data(self, adata: AnnData) -> Optional[dict]:
        """
        Returns a dict with length number of batches where each entry is a mask.

        The mask is over cell measurement columns that are present (observed)
        in each batch. Absence is defined by all 0 for that protein in that batch.
        """
        pro_exp = self.get_field(adata)
        pro_exp = pro_exp.to_numpy() if isinstance(pro_exp, pd.DataFrame) else pro_exp
        batches = adata.obs[self.batch_key].values
        batch_mask = {}
        for b_code, b in enumerate(np.unique(batches)):
            b_inds = np.where(batches.ravel() == b)[0]
            batch_sum = pro_exp[b_inds, :].sum(axis=0)
            all_zero = batch_sum == 0
            batch_mask[b_code] = ~all_zero

        if np.sum([~b[1] for b in batch_mask.items()]) > 0:
            logger.info("Found batches with missing protein expression")
            return batch_mask

        return None

    def register_field(self, adata: AnnData) -> dict:
        state_registry = super().register_field(adata)

        batch_mask = self._get_batch_mask_protein_data(adata)
        if batch_mask is not None:
            state_registry[self.PROTEIN_BATCH_MASK] = batch_mask

        return state_registry

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        transfer_state_registry = super().transfer_field(
            state_registry, adata_target, **kwargs
        )
        batch_mask = self._get_batch_mask_protein_data(adata_target)
        if batch_mask is not None:
            transfer_state_registry[self.PROTEIN_BATCH_MASK] = batch_mask

        return transfer_state_registry
