import logging
import warnings
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
import rich
from anndata import AnnData
from pandas.api.types import CategoricalDtype

from scvi.data import _constants
from scvi.data._utils import (
    _check_nonnegative_integers,
    _make_column_categorical,
    _verify_and_correct_data_format,
)

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
    def registry_key(self) -> str:
        return self._registry_key

    @property
    def attr_name(self) -> str:
        return self._attr_name


class ObsmField(BaseObsmField):
    """
    An AnnDataField for an .obsm field in the AnnData data structure.

    In addition to creating a reference to the .obsm field, stores the column
    keys for the obsm field in a more accessible .uns attribute.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obsm_key
        Key to access the field in the AnnData .obsm mapping.
    colnames_uns_key
        Key to access column names corresponding to each column of the .obsm field in
        the AnnData .uns mapping. If None, checks if the field is stored as a dataframe.
        If so, uses the dataframe's colnames. Otherwise, generates sequential column names
        (e.g. 1, 2, 3, etc.).
    is_count_data
        If True, checks if the data are counts during validation.
    correct_data_format
        If True, checks and corrects that the AnnData field is C_CONTIGUOUS and csr
        if it is dense numpy or sparse respectively.
    """

    COLUMN_NAMES_KEY = "column_names"

    def __init__(
        self,
        registry_key: str,
        obsm_key: str,
        colnames_uns_key: Optional[str] = None,
        is_count_data: bool = False,
        correct_data_format: bool = True,
    ) -> None:
        super().__init__(registry_key)
        self._attr_key = obsm_key
        self.colnames_uns_key = colnames_uns_key
        self.is_count_data = is_count_data
        self.correct_data_format = correct_data_format
        self.count_stat_key = f"n_{self.registry_key}"

    @property
    def attr_key(self) -> str:
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return False

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        if self.attr_key not in adata.obsm:
            raise KeyError(f"{self.attr_key} not found in adata.obsm.")

        obsm_data = self.get_field_data(adata)

        if self.is_count_data and not _check_nonnegative_integers(obsm_data):
            warnings.warn(
                f"adata.obsm['{self.attr_key}'] does not contain unnormalized count data. "
                "Are you sure this is what you want?"
            )

    def _setup_column_names(self, adata: AnnData) -> Union[list, np.ndarray]:
        """
        Returns a list or NumPy array of column names that will be used for the relevant .obsm data.

        If the ``colnames_uns_key`` was specified, then the columns stored in that
        field will be returned. Otherwise, if the stored data is a pandas dataframe, then
        the dataframe's colnames will be returned. In the case the stored data is a NumPy array,
        sequential column names will be generated (e.g. 1, 2, 3, etc.)
        """
        obsm_data = self.get_field_data(adata)
        if self.colnames_uns_key is None and isinstance(obsm_data, pd.DataFrame):
            logger.info(
                f"Using column names from columns of adata.obsm['{self.attr_key}']"
            )
            column_names = list(obsm_data.columns)
        elif self.colnames_uns_key is not None:
            logger.info(f"Using column names from adata.uns['{self.colnames_uns_key}']")
            column_names = adata.uns[self.colnames_uns_key]
        else:
            logger.info("Generating sequential column names")
            column_names = np.arange(obsm_data.shape[1])
        return column_names

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        if self.correct_data_format:
            _verify_and_correct_data_format(adata, self.attr_name, self.attr_key)

        column_names = self._setup_column_names(adata)

        return {self.COLUMN_NAMES_KEY: column_names}

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)
        self.validate_field(adata_target)
        source_cols = state_registry[self.COLUMN_NAMES_KEY]
        target_data = self.get_field_data(adata_target)
        if len(source_cols) != target_data.shape[1]:
            raise ValueError(
                f"Target adata.obsm['{self.attr_key}'] has {target_data.shape[1]} which does not match "
                f"the source adata.obsm['{self.attr_key}'] column count of {len(source_cols)}."
            )

        if isinstance(target_data, pd.DataFrame) and source_cols != list(
            target_data.columns
        ):
            raise ValueError(
                f"Target adata.obsm['{self.attr_key}'] column names do not match "
                f"the source adata.obsm['{self.attr_key}'] column names."
            )

        return {self.COLUMN_NAMES_KEY: state_registry[self.COLUMN_NAMES_KEY].copy()}

    def get_summary_stats(self, state_registry: dict) -> dict:
        n_obsm_cols = len(state_registry[self.COLUMN_NAMES_KEY])
        return {self.count_stat_key: n_obsm_cols}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        return None


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
            if obs_key not in adata.obs:
                raise KeyError(f"{obs_key} not found in adata.obs.")

    def _combine_obs_fields(self, adata: AnnData) -> None:
        adata.obsm[self.attr_key] = adata.obs[self.obs_keys].copy()

    @property
    def obs_keys(self) -> List[str]:
        """List of .obs keys that make up this joint field."""
        return self._obs_keys

    @property
    def attr_key(self) -> str:
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

    COLUMNS_KEY = "columns"

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key, obs_keys)

        self.count_stat_key = f"n_{self.registry_key}"

    def register_field(self, adata: AnnData) -> dict:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        return {self.COLUMNS_KEY: adata.obsm[self.attr_key].columns.to_numpy()}

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> dict:
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        n_obs_keys = len(self.obs_keys)
        return {self.count_stat_key: n_obs_keys}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        if self.is_empty:
            return None

        t = rich.table.Table(title=f"{self.registry_key} State Registry")
        t.add_column(
            "Source Location",
            justify="center",
            style="dodger_blue1",
            no_wrap=True,
            overflow="fold",
        )
        for key in state_registry[self.COLUMNS_KEY]:
            t.add_row("adata.obs['{}']".format(key))
        return t


class CategoricalJointObsField(JointObsField):
    """
    An AnnDataField for a collection of categorical .obs fields in the AnnData data structure.

    Creates an .obsm field compiling the given .obs fields. The model will reference the compiled
    data as a whole.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    MAPPINGS_KEY = "mappings"
    FIELD_KEYS_KEY = "field_keys"
    N_CATS_PER_KEY = "n_cats_per_key"

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key, obs_keys)
        self.count_stat_key = f"n_{self.registry_key}"

    def _default_mappings_dict(self) -> dict:
        return {
            self.MAPPINGS_KEY: dict(),
            self.FIELD_KEYS_KEY: [],
            self.N_CATS_PER_KEY: [],
        }

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
            categorical_dtype = (
                CategoricalDtype(categories=category_dict[key])
                if category_dict is not None
                else None
            )
            mapping = _make_column_categorical(
                obsm_df, key, key, categorical_dtype=categorical_dtype
            )
            categories[key] = mapping

        store_cats = categories if category_dict is None else category_dict

        mappings_dict = self._default_mappings_dict()
        mappings_dict[self.MAPPINGS_KEY] = store_cats
        mappings_dict[self.FIELD_KEYS_KEY] = self.obs_keys
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

    def get_summary_stats(self, _state_registry: dict) -> dict:
        n_obs_keys = len(self.obs_keys)

        return {
            self.count_stat_key: n_obs_keys,
        }

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        if self.is_empty:
            return None

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
        for key, mappings in state_registry[self.MAPPINGS_KEY].items():
            for i, mapping in enumerate(mappings):
                if i == 0:
                    t.add_row("adata.obs['{}']".format(key), str(mapping), str(i))
                else:
                    t.add_row("", str(mapping), str(i))
            t.add_row("", "")
        return t
