import logging
from typing import Optional

import rich
from anndata import AnnData

from scvi.data import _constants

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseUnsField(BaseAnnDataField):
    """
    An abstract AnnDataField for .uns attributes in the AnnData data structure.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    uns_key
        Key to access the field in the AnnData uns mapping. If None, defaults to `registry_key`.
    required
        If False, allows for `uns_key is None` and marks the field as `is_empty`.
    """

    _attr_name = _constants._ADATA_ATTRS.UNS

    def __init__(
        self, registry_key: str, uns_key: Optional[str], required: bool = True
    ) -> None:
        super().__init__()
        if required and uns_key is None:
            raise ValueError(
                "`uns_key` cannot be `None` if `required=True`. Please provide an `uns_key`."
            )
        self._registry_key = registry_key
        self._attr_key = uns_key
        self._is_empty = uns_key is None

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


class StringUnsField(BaseUnsField):
    """An AnnDataField for string .uns attributes in the AnnData data structure."""

    def validate_field(self, adata: AnnData) -> None:
        """Validate the field."""
        super().validate_field(adata)
        if self.attr_key not in adata.uns:
            raise KeyError(f"{self.attr_key} not found in adata.uns.")

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        return super().register_field(adata)

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> dict:
        """Transfer the field."""
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        """Get summary stats."""
        return {}

    def view_state_registry(self, _state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
        return None
