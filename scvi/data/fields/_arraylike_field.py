import logging
import warnings
from typing import Dict, List, Literal, Optional, Union

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
from ._mudata import MuDataWrapper

logger = logging.getLogger(__name__)


class BaseArrayLikeField(BaseAnnDataField):
    """An abstract AnnDataField for .obsm or .varm attributes in the AnnData data structure."""

    def __init__(
        self,
        registry_key: str,
    ) -> None:
        super().__init__()
        self._registry_key = registry_key
        self._attr_name = None

    @property
    def registry_key(self) -> str:  # noqa: D102
        return self._registry_key

    @property
    def attr_name(self) -> str:  # noqa: D102
        return self._attr_name


class ArrayLikeField(BaseArrayLikeField):
    """
    An AnnDataField for an .obsm or .varm field in the AnnData data structure.

    In addition to creating a reference to the .obsm or .varm field, stores the column
    keys for the obsm or varm field in a more accessible .uns attribute.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    attr_key
        Key to access the field in the AnnData .obsm or .varm mapping.
    field_type
        Type of field. Can be either "obsm" or "varm".
    colnames_uns_key
        Key to access column names corresponding to each column of the .obsm or .varm
        field in the AnnData .uns mapping. If None, checks if the field is stored as a
        dataframe. If so, uses the dataframe's colnames. Otherwise, generates sequential
        column names (e.g. 1, 2, 3, etc.).
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
        attr_key: str,
        field_type: Literal["obsm", "varm"] = None,
        colnames_uns_key: Optional[str] = None,
        is_count_data: bool = False,
        correct_data_format: bool = True,
    ) -> None:
        super().__init__(registry_key)
        if field_type == "obsm":
            self._attr_name = _constants._ADATA_ATTRS.OBSM
        elif field_type == "varm":
            self._attr_name = _constants._ADATA_ATTRS.VARM
        else:
            raise ValueError("`field_type` must be either 'obsm' or 'varm'.")

        self._attr_key = attr_key
        self.colnames_uns_key = colnames_uns_key
        self.is_count_data = is_count_data
        self.correct_data_format = correct_data_format
        self.count_stat_key = f"n_{self.registry_key}"

    @property
    def attr_key(self) -> str:  # noqa: D102
        return self._attr_key

    @property
    def is_empty(self) -> bool:  # noqa: D102
        return False

    def validate_field(self, adata: AnnData) -> None:
        """Validate the field."""
        super().validate_field(adata)
        if self.attr_key not in getattr(adata, self.attr_name):
            raise KeyError(f"{self.attr_key} not found in adata.{self.attr_name}.")

        array_data = self.get_field_data(adata)

        if self.is_count_data and not _check_nonnegative_integers(array_data):
            warnings.warn(
                f"adata.{self.attr_name}['{self.attr_key}'] does not contain "
                "unnormalized count data. Are you sure this is what you want?"
            )

    def _setup_column_names(self, adata: AnnData) -> Union[list, np.ndarray]:
        """
        Returns a list or NumPy array of column names that will be used for the relevant .obsm data.

        If the ``colnames_uns_key`` was specified, then the columns stored in that
        field will be returned. Otherwise, if the stored data is a pandas dataframe, then
        the dataframe's colnames will be returned. In the case the stored data is a NumPy array,
        sequential column names will be generated (e.g. 1, 2, 3, etc.)
        """
        array_data = self.get_field_data(adata)
        if self.colnames_uns_key is None and isinstance(array_data, pd.DataFrame):
            logger.info(
                f"Using column names from columns of adata.{self.attr_name}['{self.attr_key}']"
            )
            column_names = list(array_data.columns)
        elif self.colnames_uns_key is not None:
            logger.info(f"Using column names from adata.uns['{self.colnames_uns_key}']")
            column_names = adata.uns[self.colnames_uns_key]
        else:
            logger.info("Generating sequential column names")
            column_names = np.arange(array_data.shape[1])
        return column_names

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        super().register_field(adata)
        if self.correct_data_format:
            _verify_and_correct_data_format(adata, self.attr_name, self.attr_key)

        column_names = self._setup_column_names(adata)

        return {self.COLUMN_NAMES_KEY: column_names}

    def transfer_field(
        self, state_registry: dict, adata_target: AnnData, **kwargs
    ) -> dict:
        """Transfer the field."""
        super().transfer_field(state_registry, adata_target, **kwargs)
        self.validate_field(adata_target)
        source_cols = state_registry[self.COLUMN_NAMES_KEY]
        target_data = self.get_field_data(adata_target)
        if len(source_cols) != target_data.shape[1]:
            raise ValueError(
                f"Target adata.{self.attr_name}['{self.attr_key}'] has {target_data.shape[1]} which does not match "
                f"the source adata.{self.attr_name}['{self.attr_key}'] column count of {len(source_cols)}."
            )

        if isinstance(target_data, pd.DataFrame) and source_cols != list(
            target_data.columns
        ):
            raise ValueError(
                f"Target adata.{self.attr_name}['{self.attr_key}'] column names do not match "
                f"the source adata.{self.attr_name}['{self.attr_key}'] column names."
            )

        return {self.COLUMN_NAMES_KEY: state_registry[self.COLUMN_NAMES_KEY].copy()}

    def get_summary_stats(self, state_registry: dict) -> dict:
        """Get summary stats."""
        n_array_cols = len(state_registry[self.COLUMN_NAMES_KEY])
        return {self.count_stat_key: n_array_cols}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
        return None


class ObsmField(ArrayLikeField):
    """An AnnDataField for an .obsm field in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="obsm", **kwargs)


class VarmField(ArrayLikeField):
    """An AnnDataField for a .varm field in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="varm", **kwargs)


MuDataObsmField = MuDataWrapper(ObsmField)
MuDataVarmField = MuDataWrapper(VarmField)


class BaseJointField(BaseArrayLikeField):
    """
    An abstract AnnDataField for a collection of .obs or .var fields in the AnnData data structure.

    Creates an .obsm or .varm field containing each .obs or .var field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    attr_keys
        Sequence of keys to combine to form the obsm or varm field.
    field_type
        Type of field. Can be either 'obsm' or 'varm'.
    """

    def __init__(
        self,
        registry_key: str,
        attr_keys: Optional[List[str]],
        field_type: Literal["obsm", "varm"] = None,
    ) -> None:
        super().__init__(registry_key)
        if field_type == "obsm":
            self._source_attr_name = _constants._ADATA_ATTRS.OBS
            self._attr_name = _constants._ADATA_ATTRS.OBSM
        elif field_type == "varm":
            self._source_attr_name = _constants._ADATA_ATTRS.VAR
            self._attr_name = _constants._ADATA_ATTRS.VARM
        else:
            raise ValueError("`field_type` must be either 'obsm' or 'varm'.")
        self._attr_key = f"_scvi_{registry_key}"
        self._attr_keys = attr_keys if attr_keys is not None else []
        self._is_empty = len(self.attr_keys) == 0

    def validate_field(self, adata: AnnData) -> None:
        """Validate the field."""
        super().validate_field(adata)
        for key in self.attr_keys:
            if key not in getattr(adata, self.source_attr_name):
                raise KeyError(f"{key} not found in adata.{self.source_attr_name}.")

    def _combine_fields(self, adata: AnnData) -> None:
        """Combine the .obs or .var fields into a single .obsm or .varm field."""
        attr = getattr(adata, self.attr_name)
        source = getattr(adata, self.source_attr_name)
        attr[self.attr_key] = source[self.attr_keys].copy()

    @property
    def attr_name(self) -> str:  # noqa: D102
        return self._attr_name

    @property
    def source_attr_name(self) -> str:  # noqa: D102
        return self._source_attr_name

    @property
    def attr_keys(self) -> List[str]:
        """List of .obs or .var keys that make up this joint field."""
        return self._attr_keys

    @property
    def attr_key(self) -> str:  # noqa: D102
        return self._attr_key

    @property
    def is_empty(self) -> bool:  # noqa: D102
        return self._is_empty


class NumericalJointField(BaseJointField):
    """
    An AnnDataField for a collection of numerical .obs or .var fields in the AnnData data structure.

    Creates an .obsm or .varm field containing each .obs or .var field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    attr_keys
        Sequence of keys to combine to form the obsm or varm field.
    field_type
        Type of field. Can be either 'obsm' or 'varm'.
    """

    COLUMNS_KEY = "columns"

    def __init__(
        self,
        registry_key: str,
        attr_keys: Optional[List[str]],
        field_type: Literal["obsm", "varm"] = None,
    ) -> None:
        super().__init__(registry_key, attr_keys, field_type=field_type)

        self.count_stat_key = f"n_{self.registry_key}"

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        super().register_field(adata)
        self._combine_fields(adata)
        return {
            self.COLUMNS_KEY: getattr(adata, self.attr_name)[
                self.attr_key
            ].columns.to_numpy()
        }

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> dict:
        """Transfer the field."""
        super().transfer_field(state_registry, adata_target, **kwargs)
        return self.register_field(adata_target)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        """Get summary stats."""
        n_keys = len(self.attr_keys)
        return {self.count_stat_key: n_keys}

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
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
            t.add_row(f"adata.{self.source_attr_name}['{key}']")
        return t


class NumericalJointObsField(NumericalJointField):
    """An AnnDataField for a collection of numerical .obs fields in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="obsm", **kwargs)


class NumericalJointVarField(NumericalJointField):
    """An AnnDataField for a collection of numerical .var fields in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="varm", **kwargs)


MuDataNumericalJointObsField = MuDataWrapper(NumericalJointObsField)
MuDataNumericalJointVarField = MuDataWrapper(NumericalJointVarField)


class CategoricalJointField(BaseJointField):
    """
    An AnnDataField for a collection of categorical .obs or .var fields in the AnnData data structure.

    Creates an .obsm or .varm field compiling the given .obs or .var fields. The model
    will reference the compiled data as a whole.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    attr_keys
        Sequence of keys to combine to form the obsm or varm field.
    field_type
        Type of field. Can be either 'obsm' or 'varm'.
    """

    MAPPINGS_KEY = "mappings"
    FIELD_KEYS_KEY = "field_keys"
    N_CATS_PER_KEY = "n_cats_per_key"

    def __init__(
        self,
        registry_key: str,
        attr_keys: Optional[List[str]],
        field_type: Literal["obsm", "varm"] = None,
    ) -> None:
        super().__init__(registry_key, attr_keys, field_type=field_type)
        self.count_stat_key = f"n_{self.registry_key}"

    def _default_mappings_dict(self) -> dict:
        return {
            self.MAPPINGS_KEY: dict(),
            self.FIELD_KEYS_KEY: [],
            self.N_CATS_PER_KEY: [],
        }

    def _make_array_categorical(
        self, adata: AnnData, category_dict: Optional[Dict[str, List[str]]] = None
    ) -> dict:
        """Make the .obsm categorical."""
        if (
            self.attr_keys
            != getattr(adata, self.attr_name)[self.attr_key].columns.tolist()
        ):
            raise ValueError(
                f"Original .{self.source_attr_name} keys do not match the columns in the ",
                f"generated .{self.attr_name} field.",
            )

        categories = dict()
        df = getattr(adata, self.attr_name)[self.attr_key]
        for key in self.attr_keys:
            categorical_dtype = (
                CategoricalDtype(categories=category_dict[key])
                if category_dict is not None
                else None
            )
            mapping = _make_column_categorical(
                df, key, key, categorical_dtype=categorical_dtype
            )
            categories[key] = mapping

        store_cats = categories if category_dict is None else category_dict

        mappings_dict = self._default_mappings_dict()
        mappings_dict[self.MAPPINGS_KEY] = store_cats
        mappings_dict[self.FIELD_KEYS_KEY] = self.attr_keys
        for k in self.attr_keys:
            mappings_dict[self.N_CATS_PER_KEY].append(len(store_cats[k]))
        return mappings_dict

    def register_field(self, adata: AnnData) -> dict:
        """Register the field."""
        super().register_field(adata)
        self._combine_fields(adata)
        return self._make_array_categorical(adata)

    def transfer_field(
        self,
        state_registry: dict,
        adata_target: AnnData,
        extend_categories: bool = False,
        **kwargs,
    ) -> dict:
        """Transfer the field."""
        super().transfer_field(state_registry, adata_target, **kwargs)

        if self.is_empty:
            return

        source_cat_dict = state_registry[self.MAPPINGS_KEY].copy()
        if extend_categories:
            for key, mapping in source_cat_dict.items():
                for c in np.unique(getattr(adata_target, self.source_attr_name)[key]):
                    if c not in mapping:
                        mapping = np.concatenate([mapping, [c]])
                source_cat_dict[key] = mapping

        self.validate_field(adata_target)
        self._combine_fields(adata_target)
        return self._make_array_categorical(adata_target, category_dict=source_cat_dict)

    def get_summary_stats(self, _state_registry: dict) -> dict:
        """Get summary stats."""
        n_keys = len(self.attr_keys)

        return {
            self.count_stat_key: n_keys,
        }

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
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
                    t.add_row(
                        f"adata.{self.source_attr_name}['{key}']", str(mapping), str(i)
                    )
                else:
                    t.add_row("", str(mapping), str(i))
            t.add_row("", "")
        return t


class CategoricalJointObsField(CategoricalJointField):
    """An AnnDataField for a collection of categorical .obs fields in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="obsm", **kwargs)


class CategoricalJointVarField(CategoricalJointField):
    """An AnnDataField for a collection of categorical .var fields in the AnnData data structure."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, field_type="varm", **kwargs)


MuDataCategoricalJointObsField = MuDataWrapper(CategoricalJointObsField)
MuDataCategoricalJointVarField = MuDataWrapper(CategoricalJointVarField)
