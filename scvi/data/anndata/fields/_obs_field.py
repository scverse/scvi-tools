import logging
from typing import Optional

import numpy as np
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _make_obs_column_categorical

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseObsField(BaseAnnDataField):
    """An abstract AnnDataField for .obs attributes in the AnnData data structure."""

    _attr_name = _constants._ADATA_ATTRS.OBS

    def __init__(self, registry_key: str, obs_key: str) -> None:
        super().__init__()
        self._registry_key = registry_key
        self._attr_key = obs_key

    @property
    def registry_key(self):
        return self._registry_key

    @property
    def attr_name(self):
        return self._attr_name

    @property
    def attr_key(self):
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return False


class CategoricalObsField(BaseObsField):
    """
    An AnnDataField for categorical .obs attributes in the AnnData data structure.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_key
        Key to access the field in the AnnData obs mapping. If None, defaults to `registry_key`.
    """

    CATEGORICAL_MAPPING_KEY = "categorical_mapping"
    CM_ORIGINAL_KEY = "original_key"
    CM_MAPPING_KEY = "mapping"

    def __init__(self, registry_key: str, obs_key: Optional[str]) -> None:
        self.is_default = obs_key is None
        self._original_attr_key = obs_key or self.registry_key
        super().__init__(registry_key, f"_scvi_{self._original_attr_key}")

        self.count_stat_key = f"n_{self.registry_key}"

    def _setup_default_attr(self, adata: AnnData) -> None:
        adata.obs[self.attr_key] = np.zeros(adata.shape[0], dtype=np.int64)

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        assert self._original_attr_key in adata.obs

    def register_field(self, adata: AnnData) -> dict:
        if self.is_default:
            self._setup_default_attr(adata)

        super().register_field(adata)
        categorical_mapping = _make_obs_column_categorical(
            adata, self._original_attr_key, self.attr_key
        )
        return {self.CATEGORICAL_MAPPING_KEY: categorical_mapping}

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        extend_categories: bool = False,
        **kwargs,
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)

        if self.is_default:
            self._setup_default_attr(adata_target)

        self.validate_field(adata_target)

        mapping = state_registry[self.CATEGORICAL_MAPPING_KEY].copy()

        # extend mapping for new categories
        if extend_categories:
            for c in np.unique(self.get_field(adata_target)):
                if c not in mapping:
                    mapping = np.concatenate([mapping, [c]])
        cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
        new_mapping = _make_obs_column_categorical(
            adata_target,
            self._original_attr_key,
            self.attr_key,
            categorical_dtype=cat_dtype,
        )
        return {self.CATEGORICAL_MAPPING_KEY: new_mapping}

    def compute_summary_stats(self, state_registry: dict) -> dict:
        categorical_mapping = state_registry[self.CATEGORICAL_MAPPING_KEY]
        n_categories = len(np.unique(categorical_mapping))
        return {self.count_stat_key: n_categories}
