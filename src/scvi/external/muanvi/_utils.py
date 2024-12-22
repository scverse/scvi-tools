import warnings
from collections.abc import Iterable as IterableClass
from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import _make_column_categorical, get_anndata_attribute
from scvi.data.fields import CategoricalJointObsField
from scvi.dataloaders import ConcatDataLoader
from scvi.dataloaders._ann_dataloader import AnnDataLoader


# Class creating a new Obsm field for partially annotated layers of labels
class LabelsWithUnlabeledJointObsField(CategoricalJointObsField):
    """
    An AnnDataField for a collection of partially observed layers of labels .obs fields in the AnnData data structure.

    Creates an .obsm field compiling the given .obs fields. The model will reference the compiled
    data as a whole.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    attr_keys
        Sequence of keys to combine to form the obsm or varm field.
    unlabeled_category
        A single category to represent unlabeled cells in the data.
    """

    MAPPINGS_KEY = "mappings"
    FIELD_KEYS_KEY = "field_keys"
    N_CATS_PER_KEY = "n_cats_per_key"
    UNLABELED_CATEGORY = "unlabeled_category"

    def __init__(
        self,
        registry_key: str,
        attr_keys: list[str] | None,
        unlabeled_category: str | None,
    ) -> None:
        super().__init__(registry_key, attr_keys)
        self.count_stat_key = f"n_{self.registry_key}"
        self.unlabeled_category = unlabeled_category

    def _default_mappings_dict(self) -> dict:
        return {
            self.MAPPINGS_KEY: dict(),
            self.FIELD_KEYS_KEY: [],
            self.N_CATS_PER_KEY: [],
            self.UNLABELED_CATEGORY: [],
        }

    def _make_obsm_categorical(
        self, adata: AnnData, category_dict: dict[str, list[str]] | None = None
    ) -> dict:
        if self.attr_keys != getattr(adata, self.attr_name)[self.attr_key].columns.tolist():
            raise ValueError(
                f"Original .{self.source_attr_name} keys do not match the columns in the ",
                f"generated .{self.attr_name} field.",
            )

        categories = {}
        df = getattr(adata, self.attr_name)[self.attr_key]
        for level, key in enumerate(self.attr_keys):
            categorical_dtype = (
                CategoricalDtype(categories=category_dict[key])
                if category_dict is not None
                else None
            )
            if categorical_dtype is None:
                categorical_obs = df[key].astype("category")
            else:
                categorical_obs = df[key].astype(categorical_dtype)

            mapping = categorical_obs.cat.categories.to_numpy(copy=True)
            mapping = self._remap_unlabeled_to_final_category(mapping, level)
            cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
            mapping = _make_column_categorical(df, key, key, categorical_dtype=cat_dtype)
            categories[key] = mapping

        store_cats = categories if category_dict is None else category_dict

        mappings_dict = self._default_mappings_dict()
        mappings_dict[self.MAPPINGS_KEY] = store_cats
        mappings_dict[self.FIELD_KEYS_KEY] = self.attr_keys
        mappings_dict[self.UNLABELED_CATEGORY] = self.unlabeled_category
        for k in self.attr_keys:
            mappings_dict[self.N_CATS_PER_KEY].append(len(store_cats[k]) - 1)
        return mappings_dict

    def _remap_unlabeled_to_final_category(self, mapping: np.ndarray, level: int) -> np.ndarray:
        # Make unlabeled category the last element
        unlabeled_category = self.unlabeled_category

        # Check if the unlabeled category is in the mapping
        if unlabeled_category in mapping:
            # Find the index of the unlabeled category
            unlabeled_idx = np.where(mapping == unlabeled_category)[0][0]
            # Swap the unlabeled category with the last element
            mapping[unlabeled_idx], mapping[-1] = mapping[-1], mapping[unlabeled_idx]
        else:
            # Append the unlabeled category if it's not in the mapping
            mapping = np.append(mapping, unlabeled_category)

        return mapping

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        self._combine_fields(adata)
        state_registry = self._make_obsm_categorical(adata)
        return state_registry

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        extend_categories: bool = False,
        allow_missing_labels: bool = False,
        **kwargs,
    ) -> dict:
        """Transfer the field."""
        for level, key in enumerate(self.attr_keys):
            if (
                allow_missing_labels
                and key is not None
                and key not in list(adata_target.obs.columns)
            ):
                # Fill in original .obs attribute with unlabeled_category values.
                warnings.warn(
                    f"Missing labels key {key}. Filling in with "
                    f"unlabeled category {self.unlabeled_category}.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
                adata_target.obs[key] = self.unlabeled_category[level]

        kwargs.pop("extend_categories", None)
        transfer_state_registry = super().transfer_field(
            state_registry, adata_target, extend_categories=extend_categories, **kwargs
        )
        categories = {}
        mapping = transfer_state_registry[self.MAPPINGS_KEY]
        for level, key in enumerate(self.attr_keys):
            mapping_ = self._remap_unlabeled_to_final_category(mapping[key], level)
            categories[key] = mapping_
        store_cats = categories
        mappings_dict = self._default_mappings_dict()
        mappings_dict[self.MAPPINGS_KEY] = store_cats
        mappings_dict[self.FIELD_KEYS_KEY] = self.attr_keys
        mappings_dict[self.UNLABELED_CATEGORY] = self.unlabeled_category
        for k in self.attr_keys:
            mappings_dict[self.N_CATS_PER_KEY].append(len(store_cats[k]) - 1)
        return mappings_dict


def _get_site_code_from_category(adata_manager: AnnDataManager, category: Sequence[int | str]):
    if not isinstance(category, IterableClass) or isinstance(category, str):
        category = [category]

    site_mappings = adata_manager.get_state_registry(REGISTRY_KEYS.SITE_KEY).categorical_mapping
    site_code = []
    for cat in category:
        if cat is None:
            site_code.append(None)
            continue
        elif isinstance(cat, int) and cat < len(site_mappings):
            site_code.append(site_mappings[cat])
        elif cat not in site_mappings:
            raise ValueError(f'"{cat}" not a valid site category.')
        else:
            site_loc = np.where(site_mappings == cat)[0][0]
        site_code.append(site_loc)
    return site_code, site_mappings
