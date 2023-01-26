import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from ._arraylike_field import ObsmField
from ._layer_field import LayerField
from ._mudata import BaseMuDataWrapperClass, MuDataWrapper

logger = logging.getLogger(__name__)


class ProteinFieldMixin:
    """
    A mixin class for an protein data stored in a field of an AnnData object.

    For usage with the TotalVI model. Computes an additional mask which indicates
    where batches are missing protein data.

    Parameters
    ----------
    use_batch_mask
        If True, computes a batch mask over the data for missing protein data.
        Requires ``batch_key`` to be not None.
    batch_key
        Key corresponding to the .obs field where batch indices are stored.
        Used for computing a batch mask over the data for missing protein data.
    """

    PROTEIN_BATCH_MASK = "protein_batch_mask"

    def __init__(
        self,
        *base_field_args,
        use_batch_mask: bool = True,
        batch_field: Optional[str] = None,
        **base_field_kwargs,
    ) -> None:
        if use_batch_mask and batch_field is None:
            raise ValueError(
                "`use_batch_mask = True` requires that `batch_field is not None`. "
                "Please provide a `batch_field`."
            )
        self.use_batch_mask = use_batch_mask
        self.batch_field = batch_field
        super().__init__(
            *base_field_args,
            **base_field_kwargs,
        )

    def _get_batch_mask_protein_data(self, adata: AnnData) -> Optional[dict]:
        """
        Returns a dict with length number of batches where each entry is a mask.

        The mask is over cell measurement columns that are present (observed)
        in each batch. Absence is defined by all 0 for that protein in that batch.
        """
        pro_exp = self.get_field_data(adata)
        pro_exp = pro_exp.to_numpy() if isinstance(pro_exp, pd.DataFrame) else pro_exp
        batches = self.batch_field.get_field_data(adata)
        batch_mask = {}
        for b in np.unique(batches):
            b_inds = np.where(batches.ravel() == b)[0]
            batch_sum = pro_exp[b_inds, :].sum(axis=0)
            all_zero = batch_sum == 0
            batch_mask[str(b)] = ~all_zero

        if np.sum([~b[1] for b in batch_mask.items()]) > 0:
            logger.info("Found batches with missing protein expression")
            return batch_mask

        return None

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        state_registry = super().register_field(adata)

        if self.use_batch_mask:
            batch_mask = self._get_batch_mask_protein_data(adata)
            if batch_mask is not None:
                state_registry[self.PROTEIN_BATCH_MASK] = batch_mask

        return state_registry

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        """Transfer the field."""
        transfer_state_registry = super().transfer_field(
            state_registry, adata_target, **kwargs
        )
        batch_mask = self._get_batch_mask_protein_data(adata_target)
        if batch_mask is not None:
            transfer_state_registry[self.PROTEIN_BATCH_MASK] = batch_mask

        return transfer_state_registry


class ProteinObsmField(ProteinFieldMixin, ObsmField):
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
    use_batch_mask
        If True, computes a batch mask over the data for missing protein data.
        Requires ``batch_key`` to be not None.
    batch_key
        Key corresponding to the .obs field where batch indices are stored.
        Used for computing a batch mask over the data for missing protein data.
    colnames_uns_key
        Key to access column names corresponding to each column of the .obsm field in
        the AnnData .uns mapping. Only used when .obsm data is a np.ndarray, not a pd.DataFrame.
    is_count_data
        If True, checks if the data are counts during validation.
    correct_data_format
        If True, checks and corrects that the AnnData field is C_CONTIGUOUS and csr
        if it is dense numpy or sparse respectively.
    """


class ProteinLayerField(ProteinFieldMixin, LayerField):
    """
    An AnnDataField for an protein data stored in a layer field of an AnnData object.

    For usage with the TotalVI model. Computes an additional mask which indicates
    where batches are missing protein data.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    layer
        Key to access the field in the AnnData layers mapping. If None, uses the data in .X.
    use_batch_mask
        If True, computes a batch mask over the data for missing protein data.
        Requires ``batch_key`` to be not None.
    batch_key
        Key corresponding to the .obs field where batch indices are stored.
        Used for computing a batch mask over the data for missing protein data.
    is_count_data
        If True, checks if the data are counts during validation.
    correct_data_format
        If True, checks and corrects that the AnnData field is C_CONTIGUOUS and csr
        if it is dense numpy or sparse respectively.
    """


def copy_over_batch_attr(self, mdata: MuData):
    """Copy over batch attributes from the original MuData object."""
    # Assign self.batch_field if not yet assigned to MuDataWrapped field.
    # Then, reassign self.adata_field.batch_field to the batch AnnDataField.
    if isinstance(self.adata_field.batch_field, BaseMuDataWrapperClass):
        self.batch_field = self.adata_field.batch_field
        self.adata_field.batch_field = self.batch_field.adata_field

    # Copy over batch data to the protein modality.
    batch_data = self.batch_field.get_field_data(mdata)
    bdata = self.get_modality(mdata)
    bdata_attr = getattr(bdata, self.batch_field.attr_name)
    bdata_attr[self.batch_field.attr_key] = batch_data


MuDataProteinLayerField = MuDataWrapper(
    ProteinLayerField, preregister_fn=copy_over_batch_attr
)
