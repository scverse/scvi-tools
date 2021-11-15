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
        self._columns_key = f"{self.registry_key}_keys"

    def register_field(self, adata: AnnData) -> None:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        adata.uns[_constants._SETUP_DICT_KEY][self._columns_key] = adata.obsm[
            self.attr_key
        ].columns.to_numpy()

    def transfer_field(
        self,
        setup_dict: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> None:
        super().transfer_field(setup_dict, adata_target, **kwargs)
        self.register_field(adata_target)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        super().compute_summary_stats(adata)
        n_continuous_covariates = len(self.obs_keys)
        stat_name = f"n_{self.registry_key}"
        return {stat_name: n_continuous_covariates}


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

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key, obs_keys)
        self._mappings_key = registry_key

    def _initialize_mappings_dict(self, adata: AnnData) -> None:
        adata.uns[_constants._SETUP_DICT_KEY][self._mappings_key] = dict(
            mappings=dict(), keys=[]
        )

    def _make_obsm_categorical(
        self, adata: AnnData, category_dict: Optional[Dict[str, List[str]]] = None
    ) -> None:
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
        mappings_dict[_constants._JO_CM_MAPPINGS_KEY] = store_cats
        mappings_dict[_constants._JO_CM_KEYS_KEY] = self.obs_keys
        n_cats_per_key = []
        for k in self.obs_keys:
            n_cats_per_key.append(len(store_cats[k]))
        mappings_dict[_constants._JO_CM_N_CATS_PER_KEY] = n_cats_per_key

    def register_field(self, adata: AnnData) -> None:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        self._initialize_mappings_dict(adata)
        self._make_obsm_categorical(adata)

    def transfer_field(
        self,
        setup_dict: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> None:
        super().transfer_field(setup_dict, adata_target, **kwargs)

        if self.is_empty:
            return

        extend_categories = getattr(kwargs, "extend_categories", False)

        source_cat_dict = setup_dict[self._mappings_key][
            _constants._JO_CM_MAPPINGS_KEY
        ].copy()
        if extend_categories:
            for key, mapping in source_cat_dict.items():
                for c in np.unique(adata_target.obs[key]):
                    if c not in mapping:
                        mapping = np.concatenate([mapping, [c]])
                source_cat_dict[key] = mapping

        self.validate_field(adata_target)
        self._combine_obs_fields(adata_target)
        self._initialize_mappings_dict(adata_target)
        self._make_obsm_categorical(adata_target, category_dict=source_cat_dict)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        super().compute_summary_stats(adata)
        n_categorical_covariates = len(self.obs_keys)
        stat_name = f"n_{self.registry_key}"

        return {
            stat_name: n_categorical_covariates,
        }
