import logging
from typing import Optional

import numpy as np
import rich
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data import _constants
from scvi.data._utils import _make_column_categorical, get_anndata_attribute

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseObsField(BaseAnnDataField):
    """
    An abstract AnnDataField for .obs attributes in the AnnData data structure.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_key
        Key to access the field in the AnnData obs mapping. If None, defaults to `registry_key`.
    required
        If False, allows for `obs_key is None` and marks the field as `is_empty`.
    """

    _attr_name = _constants._ADATA_ATTRS.OBS

    def __init__(
        self, registry_key: str, obs_key: Optional[str], required: bool = True
    ) -> None:
        super().__init__()
        if required and obs_key is None:
            raise ValueError(
                "`obs_key` cannot be `None` if `required=True`. Please provide an `obs_key`."
            )
        self._registry_key = registry_key
        self._attr_key = obs_key
        self._is_empty = obs_key is None

    @property
    def registry_key(self) -> str:
        return self._registry_key

    @property
    def attr_name(self) -> str:
        return self._attr_name

    @property
    def attr_key(self) -> str:
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return self._is_empty


class NumericalObsField(BaseObsField):
    """
    An AnnDataField for numerical .obs attributes in the AnnData data structure.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_key
        Key to access the field in the AnnData obs mapping. If None, defaults to `registry_key`.
    """

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        if self.attr_key not in adata.obs:
            raise KeyError(f"{self.attr_key} not found in adata.obs.")

    def register_field(self, adata: AnnData) -> dict:
        return super().register_field(adata)

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        return {}

    def view_state_registry(self, _state_registry: dict) -> Optional[rich.table.Table]:
        return None


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
    ORIGINAL_ATTR_KEY = "original_key"

    def __init__(self, registry_key: str, obs_key: Optional[str]) -> None:
        self.is_default = obs_key is None
        self._original_attr_key = obs_key or registry_key
        super().__init__(registry_key, f"_scvi_{registry_key}")

        self.count_stat_key = f"n_{self.registry_key}"

    def _setup_default_attr(self, adata: AnnData) -> None:
        self._original_attr_key = self.attr_key
        adata.obs[self.attr_key] = np.zeros(adata.shape[0], dtype=np.int64)

    def _get_original_column(self, adata: AnnData) -> np.ndarray:
        return get_anndata_attribute(adata, self.attr_name, self._original_attr_key)

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        if self._original_attr_key not in adata.obs:
            raise KeyError(f"{self._original_attr_key} not found in adata.obs.")

    def register_field(self, adata: AnnData) -> dict:
        if self.is_default:
            self._setup_default_attr(adata)

        super().register_field(adata)
        categorical_mapping = _make_column_categorical(
            adata.obs,
            self._original_attr_key,
            self.attr_key,
        )
        return {
            self.CATEGORICAL_MAPPING_KEY: categorical_mapping,
            self.ORIGINAL_ATTR_KEY: self._original_attr_key,
        }

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
        new_mapping = _make_column_categorical(
            adata_target.obs,
            self._original_attr_key,
            self.attr_key,
            categorical_dtype=cat_dtype,
        )
        return {
            self.CATEGORICAL_MAPPING_KEY: new_mapping,
            self.ORIGINAL_ATTR_KEY: self._original_attr_key,
        }

    def get_summary_stats(self, state_registry: dict) -> dict:
        categorical_mapping = state_registry[self.CATEGORICAL_MAPPING_KEY]
        n_categories = len(np.unique(categorical_mapping))
        return {self.count_stat_key: n_categories}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        source_key = state_registry[self.ORIGINAL_ATTR_KEY]
        mapping = state_registry[self.CATEGORICAL_MAPPING_KEY]
        t = rich.table.Table(title=f"{self.registry_key} State Registry")
        t.add_column(
            "Source Location",
            justify="center",
            style="dodger_blue1",
            no_wrap=True,
            overflow="fold",
        )
        t.add_column(
            "Categories", justify="center", style="green", no_wrap=True, overflow="fold"
        )
        t.add_column(
            "scvi-tools Encoding",
            justify="center",
            style="dark_violet",
            no_wrap=True,
            overflow="fold",
        )
        for i, cat in enumerate(mapping):
            if i == 0:
                t.add_row("adata.obs['{}']".format(source_key), str(cat), str(i))
            else:
                t.add_row("", str(cat), str(i))
        return t
