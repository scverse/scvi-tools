from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from anndata import AnnData

from . import _constants
from ._utils import _assert_key_in_obs, _make_obs_column_categorical


class BaseAnnDataField(ABC):
    def __init__(self) -> None:
        super().__init__()
        self.registered = False

    @abstractmethod
    def register_field(self, adata: AnnData) -> None:
        self.validate(adata)
        self.registered = True

    @abstractmethod
    def validate_field(self, adata: AnnData) -> None:
        pass

    @abstractmethod
    def data_registry_mapping(self) -> dict:
        pass

    def compute_summary_stats(self, adata: AnnData) -> dict:
        assert self.registered
        return dict()


class BaseObsField(BaseAnnDataField):
    attr_name = _constants._ADATA_ATTRS.OBS

    def __init__(self, scvi_key: str, attr_key: str) -> None:
        super().__init__()
        self.scvi_key = scvi_key
        self.attr_key = attr_key

    def validate_field(self, adata: AnnData) -> None:
        super().validate(adata)
        _assert_key_in_obs(adata, self.attr_key)


class CategoricalObsField(BaseObsField):
    def __init__(self, scvi_key: str, attr_key: Optional[str]) -> None:
        super().__init__(scvi_key, attr_key)
        self.is_default = attr_key is None
        self.attr_key = attr_key or self.scvi_key
        self.category_code_key = f"{self.attr_key}_scvi"

    def _setup_default_attr(self, adata: AnnData) -> None:
        adata.obs[self.attr_key] = np.zeros(adata.shape[0], dtype=np.int64)

    def register_field(self, adata: AnnData) -> None:
        if self.is_default:
            self._setup_default_attr(adata)

        super().register_field(adata)
        _make_obs_column_categorical(
            adata, self.attr_key, alternate_column_key=self.category_code_key
        )

    def data_registry_mapping(self) -> dict:
        super().data_registry_mapping()
        return {
            self.scvi_key: {
                _constants._DR_ATTR_NAME: self.attr_name,
                _constants._DR_ATTR_KEY: self.category_code_key,
            }
        }

    def compute_summary_stats(self, adata: AnnData) -> dict:
        super().get_summary_stats(adata)
        categorical_mappings = self.adata.uns[_constants._SETUP_DICT_KEY][
            _constants._CATEGORICAL_MAPPINGS_KEY
        ]
        n_categories = len(
            np.unique(categorical_mappings[self.scvi_key][_constants._CM_MAPPING_KEY])
        )
        stat_name = f"n_{self.scvi_key}"
        return {stat_name: n_categories}
