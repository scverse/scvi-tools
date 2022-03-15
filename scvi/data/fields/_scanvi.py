from typing import Optional, Union

import numpy as np
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data._utils import _make_column_categorical

from ._obs_field import CategoricalObsField


class LabelsWithUnlabeledObsField(CategoricalObsField):
    """
    An AnnDataField for labels which include explicitly unlabeled cells.

    Remaps the unlabeled category to the final index if present in labels.
    The unlabeled category is a specific category name specified by the user.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_key
        Key to access the field in the AnnData obs mapping. If None, defaults to `registry_key`.
    unlabeled_category
        Value assigned to unlabeled cells.
    """

    UNLABELED_CATEGORY = "unlabeled_category"

    def __init__(
        self,
        registry_key: str,
        obs_key: Optional[str],
        unlabeled_category: Union[str, int, float],
    ) -> None:
        super().__init__(registry_key, obs_key)
        self._unlabeled_category = unlabeled_category

    def _remap_unlabeled_to_final_category(
        self, adata: AnnData, mapping: np.ndarray
    ) -> dict:
        labels = self._get_original_column(adata)

        if self._unlabeled_category in labels:
            unlabeled_idx = np.where(mapping == self._unlabeled_category)
            unlabeled_idx = unlabeled_idx[0][0]
            # move unlabeled category to be the last position
            mapping[unlabeled_idx], mapping[-1] = mapping[-1], mapping[unlabeled_idx]
        # could be in mapping in transfer case
        elif self._unlabeled_category not in mapping:
            # just put as last category
            mapping = np.asarray(list(mapping) + [self._unlabeled_category])

        cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
        # rerun setup for the batch column
        adata = self._maybe_get_modality(adata)
        mapping = _make_column_categorical(
            adata.obs,
            self._original_attr_key,
            self.attr_key,
            categorical_dtype=cat_dtype,
        )

        return {
            self.CATEGORICAL_MAPPING_KEY: mapping,
            self.ORIGINAL_ATTR_KEY: self._original_attr_key,
            self.UNLABELED_CATEGORY: self._unlabeled_category,
        }

    def register_field(self, adata: AnnData) -> dict:
        if self.is_default:
            self._setup_default_attr(adata)

        state_registry = super().register_field(adata)
        mapping = state_registry[self.CATEGORICAL_MAPPING_KEY]
        return self._remap_unlabeled_to_final_category(adata, mapping)

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        extend_categories: bool = False,
        **kwargs,
    ) -> dict:
        transfer_state_registry = super().transfer_field(
            state_registry, adata_target, extend_categories=extend_categories, **kwargs
        )
        mapping = transfer_state_registry[self.CATEGORICAL_MAPPING_KEY]
        return self._remap_unlabeled_to_final_category(adata_target, mapping)
