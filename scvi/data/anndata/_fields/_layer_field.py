import warnings
from typing import Optional

from anndata import AnnData

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _check_nonnegative_integers

from ._base_field import BaseAnnDataField


class LayerField(BaseAnnDataField):
    def __init__(self, scvi_key: str, layer: Optional[str]) -> None:
        super().__init__()
        self._scvi_key = scvi_key
        self._attr_name = (
            _constants._ADATA_ATTRS.X
            if layer is None
            else _constants._ADATA_ATTRS.LAYERS
        )
        self._attr_key = layer

    @property
    def scvi_key(self):
        return self._scvi_key

    @property
    def attr_name(self):
        return self._attr_name

    @property
    def attr_key(self):
        return self._attr_key

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        x = self.get_field(adata)

        if _check_nonnegative_integers(x) is False:
            logger_data_loc = (
                "adata.X" if self.attr_key is None else f"adata.layers[{self.attr_key}]"
            )
            warnings.warn(
                f"{logger_data_loc} does not contain unnormalized count data. "
                "Are you sure this is what you want?"
            )

    def register_field(self, adata: AnnData) -> None:
        super().register_field(adata)

    def transfer_field(self, adata_source: AnnData, adata_target: AnnData) -> None:
        super().transfer_field(adata_source, adata_target)
        summary_stats = adata_source.uns[_constants._SETUP_DICT_KEY][
            _constants._SUMMARY_STATS_KEY
        ]
        target_n_vars = adata_target.n_vars
        if target_n_vars != summary_stats["n_vars"]:
            raise ValueError(
                "Number of vars in adata_target not the same as source. "
                + "Expected: {} Received: {}".format(
                    target_n_vars, summary_stats["n_vars"]
                )
            )

        self.register_field(adata_target)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        return dict(n_cells=adata.n_obs, n_vars=adata.n_vars)
