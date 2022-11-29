import logging
import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd
import rich
from anndata import AnnData

from scvi.data import _constants
from scvi.data._utils import (
    _check_nonnegative_integers,
    _verify_and_correct_data_format,
)

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseVarmField(BaseAnnDataField):
    """An abstract AnnDataField for .varm attributes in the AnnData data structure."""

    _attr_name = _constants._ADATA_ATTRS.VARM

    def __init__(
        self,
        registry_key: str,
    ) -> None:
        super().__init__()
        self._registry_key = registry_key

    @property
    def registry_key(self) -> str:  # noqa: D102
        return self._registry_key

    @property
    def attr_name(self) -> str:  # noqa: D102
        return self._attr_name


class VarmField(BaseVarmField):
    """
    An AnnDataField for an .varm field in the AnnData data structure.
    In addition to creating a reference to the .varm field, stores the column
    keys for the varm field in a more accessible .uns attribute.
    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    varm_key
        Key to access the field in the AnnData .varm mapping.
    colnames_uns_key
        Key to access column names corresponding to each column of the .varm field in
        the AnnData .uns mapping. If None, checks if the field is stored as a dataframe.
        If so, uses the dataframe's colnames. Otherwise, generates sequential column names
        (e.g. 1, 2, 3, etc.).
    is_count_data
        If True, checks if the data are counts during validation.
    correct_data_format
        If True, checks and corrects that the AnnData field is C_CONTIGUOUS and csr
        if it is dense numpy or sparse respectively.
    """

    COLUMN_NAMES_KEY = "column_names"

    def __init__(
        self,
        registry_key: str,
        varm_key: str,
        colnames_uns_key: Optional[str] = None,
        is_count_data: bool = False,
        correct_data_format: bool = True,
    ) -> None:
        super().__init__(registry_key)
        self._attr_key = varm_key
        self.colnames_uns_key = colnames_uns_key
        self.is_count_data = is_count_data
        self.correct_data_format = correct_data_format
        self.count_stat_key = f"n_{self.registry_key}"

    @property
    def attr_key(self) -> str:  # noqa: D102
        return self._attr_key

    @property
    def is_empty(self) -> bool:  # noqa: D102
        return False

    def validate_field(self, adata: AnnData) -> None:
        """Validate the field."""
        super().validate_field(adata)
        if self.attr_key not in adata.varm:
            raise KeyError(f"{self.attr_key} not found in adata.varm.")

        varm_data = self.get_field_data(adata)

        if self.is_count_data and not _check_nonnegative_integers(varm_data):
            warnings.warn(
                f"adata.varm['{self.attr_key}'] does not contain unnormalized count data. "
                "Are you sure this is what you want?"
            )

    def _setup_column_names(self, adata: AnnData) -> Union[list, np.ndarray]:
        """
        Returns a list or NumPy array of column names that will be used for the relevant .varm data.
        If the ``colnames_uns_key`` was specified, then the columns stored in that
        field will be returned. Otherwise, if the stored data is a pandas dataframe, then
        the dataframe's colnames will be returned. In the case the stored data is a NumPy array,
        sequential column names will be generated (e.g. 1, 2, 3, etc.)
        """
        varm_data = self.get_field_data(adata)
        if self.colnames_uns_key is None and isinstance(varm_data, pd.DataFrame):
            logger.info(
                f"Using column names from columns of adata.varm['{self.attr_key}']"
            )
            column_names = list(varm_data.columns)
        elif self.colnames_uns_key is not None:
            logger.info(f"Using column names from adata.uns['{self.colnames_uns_key}']")
            column_names = adata.uns[self.colnames_uns_key]
        else:
            logger.info("Generating sequential column names")
            column_names = np.arange(varm_data.shape[1])
        return column_names

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        super().register_field(adata)
        if self.correct_data_format:
            _verify_and_correct_data_format(adata, self.attr_name, self.attr_key)

        column_names = self._setup_column_names(adata)

        return {self.COLUMN_NAMES_KEY: column_names}

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        """Transfer the field."""
        super().transfer_field(state_registry, adata_target, **kwargs)
        self.validate_field(adata_target)
        source_cols = state_registry[self.COLUMN_NAMES_KEY]
        target_data = self.get_field_data(adata_target)
        if len(source_cols) != target_data.shape[1]:
            raise ValueError(
                f"Target adata.varm['{self.attr_key}'] has {target_data.shape[1]} which does not match "
                f"the source adata.varm['{self.attr_key}'] column count of {len(source_cols)}."
            )

        if isinstance(target_data, pd.DataFrame) and source_cols != list(
            target_data.columns
        ):
            raise ValueError(
                f"Target adata.varm['{self.attr_key}'] column names do not match "
                f"the source adata.varm['{self.attr_key}'] column names."
            )

        return {self.COLUMN_NAMES_KEY: state_registry[self.COLUMN_NAMES_KEY].copy()}

    def get_summary_stats(self, state_registry: dict) -> dict:
        """Get summary stats."""
        n_varm_cols = len(state_registry[self.COLUMN_NAMES_KEY])
        return {self.count_stat_key: n_varm_cols}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
        return None
