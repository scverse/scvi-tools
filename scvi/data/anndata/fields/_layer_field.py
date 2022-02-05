import warnings
from typing import Optional

import rich
from anndata import AnnData

from scvi.data._utils import _check_nonnegative_integers
from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _verify_and_correct_data_format

from ._base_field import BaseAnnDataField


class LayerField(BaseAnnDataField):
    """
    An AnnDataField for layer or X attributes in the AnnData data structure.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    layer
        Key to access the field in the AnnData layers mapping. If None, uses the data in .X.
    is_count_data
        If True, checks if the data are counts during validation.
    correct_data_format
        If True, checks and corrects that the AnnData field is C_CONTIGUOUS and csr
        if it is dense numpy or sparse respectively.
    """

    N_VARS_KEY = "n_vars"
    N_CELLS_KEY = "n_cells"

    def __init__(
        self,
        registry_key: str,
        layer: Optional[str],
        is_count_data: bool = True,
        correct_data_format: bool = True,
    ) -> None:
        super().__init__()
        self._registry_key = registry_key
        self._attr_name = (
            _constants._ADATA_ATTRS.X
            if layer is None
            else _constants._ADATA_ATTRS.LAYERS
        )
        self._attr_key = layer
        self.is_count_data = is_count_data
        self.correct_data_format = correct_data_format

    @property
    def registry_key(self) -> str:
        return self._registry_key

    @property
    def attr_name(self) -> str:
        return self._attr_name

    @property
    def attr_key(self) -> Optional[str]:
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return False

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        x = self.get_field_data(adata)

        if self.is_count_data and not _check_nonnegative_integers(x):
            logger_data_loc = (
                "adata.X" if self.attr_key is None else f"adata.layers[{self.attr_key}]"
            )
            warnings.warn(
                f"{logger_data_loc} does not contain unnormalized count data. "
                "Are you sure this is what you want?"
            )

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        if self.correct_data_format:
            _verify_and_correct_data_format(adata, self.attr_name, self.attr_key)
        return {self.N_CELLS_KEY: adata.n_obs, self.N_VARS_KEY: adata.n_vars}

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)
        n_vars = state_registry[self.N_VARS_KEY]
        target_n_vars = adata_target.n_vars
        if target_n_vars != n_vars:
            raise ValueError(
                "Number of vars in adata_target not the same as source. "
                + "Expected: {} Received: {}".format(target_n_vars, n_vars)
            )

        return self.register_field(adata_target)

    def get_summary_stats(self, state_registry: dict) -> dict:
        return state_registry.copy()

    def view_state_registry(self, _state_registry: dict) -> Optional[rich.table.Table]:
        return None
