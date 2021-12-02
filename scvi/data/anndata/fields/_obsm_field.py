import logging
from typing import Dict, List, Optional

import numpy as np
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data.anndata import _constants

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseObsmField(BaseAnnDataField):
    """An abstract AnnDataField for .obsm attributes in the AnnData data structure."""

    _attr_name = _constants._ADATA_ATTRS.OBSM

    def __init__(
        self,
        registry_key: str,
    ) -> None:
        super().__init__()
        self._registry_key = registry_key

    @property
    def registry_key(self):
        return self._registry_key

    @property
    def attr_name(self):
        return self._attr_name


class JointObsField(BaseObsmField):
    """
    An abstract AnnDataField for a collection of .obs fields in the AnnData data structure.

    Creates an .obsm field containing each .obs field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key)
        self._attr_key = f"_scvi_{registry_key}"
        self._obs_keys = obs_keys if obs_keys is not None else []
        self._is_empty = len(self.obs_keys) == 0

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        for obs_key in self._obs_keys:
            assert obs_key in adata.obs, f"{obs_key} not found in adata.obs."

    def _combine_obs_fields(self, adata: AnnData) -> None:
        adata.obsm[self.attr_key] = adata.obs[self.obs_keys].copy()

    @property
    def obs_keys(self):
        return self._obs_keys

    @property
    def attr_key(self):
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return self._is_empty


class NumericalJointObsField(JointObsField):
    """
    An AnnDataField for a collection of numerical .obs fields in the AnnData data structure.

    Creates an .obsm field containing each .obs field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key, obs_keys)
        self.columns_key = f"{self.registry_key}_keys"

        self.count_stat_key = f"n_{self.registry_key}"

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        return {self._columns_key: adata.obsm[self.attr_key].columns.to_numpy()}

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def compute_summary_stats(self, state_registry: dict) -> dict:
        n_continuous_covariates = len(state_registry[self.columns_key].shape[0])
        return {self.count_stat_key: n_continuous_covariates}


class CategoricalJointObsField(JointObsField):
    """
    An AnnDataField for a collection of categorical .obs fields in the AnnData data structure.

    Creates an .obsm field containing each .obs field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    MAPPINGS_KEY = "mappings"
    KEYS_KEY = "keys"
    N_CATS_PER_KEY = "n_cats_per_key"

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key, obs_keys)
        self.count_stat_key = f"n_{self.registry_key}"

    def _default_mappings_dict(self) -> dict:
        return {self.MAPPINGS_KEY: dict(), self.KEYS_KEY: [], self.N_CATS_PER_KEY: []}

    def _make_obsm_categorical(
        self, adata: AnnData, category_dict: Optional[Dict[str, List[str]]] = None
    ) -> dict:
        if self.obs_keys != adata.obsm[self.attr_key].columns.tolist():
            raise ValueError(
                "Original .obs keys do not match the columns in the generated .obsm field."
            )

        categories = dict()
        obsm_df = adata.obsm[self.attr_key]
        for key in self.obs_keys:
            if category_dict is None:
                categorical_obs = obsm_df[key].astype("category")
                mapping = categorical_obs.cat.categories.to_numpy(copy=True)
                categories[key] = mapping
            else:
                possible_cats = category_dict[key]
                categorical_obs = obsm_df[key].astype(
                    CategoricalDtype(categories=possible_cats)
                )
            obsm_df[key] = categorical_obs.cat.codes

        store_cats = categories if category_dict is None else category_dict

        mappings_dict = adata.uns[_constants._SETUP_DICT_KEY][self._mappings_key]
        mappings_dict = self._default_mappings_dict()
        mappings_dict[self.MAPPINGS_KEY] = store_cats
        mappings_dict[self.KEYS_KEY] = self.obs_keys
        for k in self.obs_keys:
            mappings_dict[self.N_CATS_PER_KEY].append(len(store_cats[k]))
        return mappings_dict

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        return self._make_obsm_categorical(adata)

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        extend_categories: bool = False,
        **kwargs,
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)

        if self.is_empty:
            return

        source_cat_dict = state_registry[self.MAPPINGS_KEY].copy()
        if extend_categories:
            for key, mapping in source_cat_dict.items():
                for c in np.unique(adata_target.obs[key]):
                    if c not in mapping:
                        mapping = np.concatenate([mapping, [c]])
                source_cat_dict[key] = mapping

        self.validate_field(adata_target)
        self._combine_obs_fields(adata_target)
        return self._make_obsm_categorical(adata_target, category_dict=source_cat_dict)

    def compute_summary_stats(self, state_registry: dict) -> dict:
        n_categorical_covariates = len(state_registry[self.KEYS_KEY])

        return {
            self.count_stat_key: n_categorical_covariates,
        }
