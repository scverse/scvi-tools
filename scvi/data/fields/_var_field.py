import logging
from typing import Optional

import numpy as np
import rich
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data import _constants
from scvi.data._utils import _make_column_categorical, get_anndata_attribute

from ._base_field import BaseAnnDataField
from ._mudata import MuDataWrapper

logger = logging.getLogger(__name__)


class BaseVarField(BaseAnnDataField):
    """
    An abstract AnnDataField for .var attributes in the AnnData data structure.
    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    var_key
        Key to access the field in the AnnData var mapping. If None, defaults to `registry_key`.
    required
        If False, allows for `var_key is None` and marks the field as `is_empty`.
    """

    _attr_name = _constants._ADATA_ATTRS.VAR

    def __init__(
        self, registry_key: str, var_key: Optional[str], required: bool = True
    ) -> None:
        super().__init__()
        if required and var_key is None:
            raise ValueError(
                "`var_key` cannot be `None` if `required=True`. Please provide an `var_key`."
            )
        self._registry_key = registry_key
        self._attr_key = var_key
        self._is_empty = var_key is None

    @property
    def registry_key(self) -> str:  # noqa: D102
        return self._registry_key

    @property
    def attr_name(self) -> str:  # noqa: D102
        return self._attr_name

    @property
    def attr_key(self) -> str:  # noqa: D102
        return self._attr_key

    @property
    def is_empty(self) -> bool:  # noqa: D102
        return self._is_empty


class NumericalVarField(BaseVarField):
    """
    An AnnDataField for numerical .var attributes in the AnnData data structure.
    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    var_key
        Key to access the field in the AnnData var mapping. If None, defaults to `registry_key`.
    """

    def validate_field(self, adata: AnnData) -> None:
        """Validate field."""
        super().validate_field(adata)
        if self.attr_key not in adata.var:
            raise KeyError(f"{self.attr_key} not found in adata.var.")

    def register_field(self, adata: AnnData) -> dict:
        """Register field."""
        return super().register_field(adata)

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        """Transfer field from registry to target AnnData."""
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        """Get summary stats."""
        return {}

    def view_state_registry(self, _state_registry: dict) -> Optional[rich.table.Table]:
        """View state registry."""
        return None
