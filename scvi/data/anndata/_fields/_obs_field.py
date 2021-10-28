from typing import Optional

import numpy as np
from anndata import AnnData

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _make_obs_column_categorical

from ._base_field import BaseAnnDataField


class BaseObsField(BaseAnnDataField):
    _attr_name = _constants._ADATA_ATTRS.OBS

    def __init__(self, scvi_key: str, obs_key: str) -> None:
        super().__init__()
        self._scvi_key = scvi_key
        self._attr_key = obs_key

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
        self.get_field(adata)


class CategoricalObsField(BaseObsField):
    def __init__(self, scvi_key: str, obs_key: Optional[str]) -> None:
        super().__init__(scvi_key, obs_key)
        self.is_default = obs_key is None
        self._attr_key = obs_key or self.scvi_key
        self.category_code_key = f"_scvi_{self.attr_key}"

    def _setup_default_attr(self, adata: AnnData) -> None:
        adata.obs[self.attr_key] = np.zeros(adata.shape[0], dtype=np.int64)

    def register_field(self, adata: AnnData) -> None:
        if self.is_default:
            self._setup_default_attr(adata)

        super().register_field(adata)
        _make_obs_column_categorical(
            adata, self.attr_key, alternate_column_key=self.category_code_key
        )

    def compute_summary_stats(self, adata: AnnData) -> dict:
        super().compute_summary_stats(adata)
        categorical_mappings = adata.uns[_constants._SETUP_DICT_KEY][
            _constants._CATEGORICAL_MAPPINGS_KEY
        ]
        n_categories = len(
            np.unique(
                categorical_mappings[self.category_code_key][_constants._CM_MAPPING_KEY]
            )
        )
        stat_name = f"n_{self.scvi_key}"
        return {stat_name: n_categories}
