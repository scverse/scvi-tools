import logging
from typing import Optional

import numpy as np
from anndata import AnnData
from mudata import MuData
from pandas.api.types import CategoricalDtype

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import (
    _make_obs_column_categorical,
    get_anndata_attribute,
    parse_attr_key,
)

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseObsField(BaseAnnDataField):
    """An abstract AnnDataField for .obs attributes in the AnnData data structure."""

    _attr_name = _constants._ADATA_ATTRS.OBS

    def __init__(self, registry_key: str, obs_key: str) -> None:
        super().__init__()
        self._registry_key = registry_key
        self._mod_key, self._attr_key = parse_attr_key(obs_key)

    @property
    def registry_key(self):
        return self._registry_key

    @property
    def mod_key(self):
        return self._mod_key

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
        super().__init__(registry_key, obs_key or registry_key)
        self.is_default = obs_key is None
        self._original_attr_key = self.attr_key
        self._original_mod_key = self.mod_key

        self._mod_key = None
        self._attr_key = f"_scvi_{self.attr_key}"

        self.count_stat_key = f"n_{self.registry_key}"

    def _setup_default_attr(self, adata: AnnData) -> None:
        adata.obs[self.attr_key] = np.zeros(adata.shape[0], dtype=np.int64)

    def _get_original_column(self, adata: AnnData) -> np.ndarray:
        return get_anndata_attribute(
            adata, self._original_mod_key, self.attr_name, self._original_attr_key
        )

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        if self._original_mod_key is not None:
            if not isinstance(adata, MuData):
                raise AssertionError("Cannot use mod_key with AnnData.")
            if self._original_mod_key not in adata.mod:
                raise AssertionError(
                    f"{self._original_mod_key} is not a valid modality."
                )
            adata = adata[self._original_mod_key]
        if self._original_attr_key not in adata.obs:
            raise AssertionError(f"{self._original_attr_key} not found in .obs.")

    def register_field(self, adata: AnnData) -> dict:
        if self.is_default:
            self._setup_default_attr(adata)

        super().register_field(adata)
        categorical_mapping = _make_obs_column_categorical(
            adata, self._original_attr_key, self.attr_key, return_mapping=True
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
        for c in np.unique(self._get_original_column(adata_target)):
            if c not in mapping:
                if extend_categories:
                    mapping = np.concatenate([mapping, [c]])
                else:
                    raise ValueError(
                        f"Category {c} not found in source registry. "
                        f"Cannot transfer setup without `extend_categories = True`."
                    )
        cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
        new_mapping = _make_obs_column_categorical(
            adata_target,
            self._original_attr_key,
            self.attr_key,
            categorical_dtype=cat_dtype,
            return_mapping=True,
        )
        return {self.CATEGORICAL_MAPPING_KEY: new_mapping}

    def get_summary_stats(self, state_registry: dict) -> dict:
        categorical_mapping = state_registry[self.CATEGORICAL_MAPPING_KEY]
        n_categories = len(np.unique(categorical_mapping))
        return {self.count_stat_key: n_categories}
