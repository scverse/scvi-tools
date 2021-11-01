from abc import ABC, abstractmethod

import numpy as np
from anndata import AnnData

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _get_field


class BaseAnnDataField(ABC):
    def __init__(self) -> None:
        super().__init__()

    @property
    @abstractmethod
    def scvi_key(self):
        pass

    @property
    @abstractmethod
    def attr_name(self):
        pass

    @property
    @abstractmethod
    def attr_key(self):
        pass

    @abstractmethod
    def validate_field(self, adata: AnnData) -> None:
        pass

    @abstractmethod
    def register_field(self, adata: AnnData) -> None:
        self.validate_field(adata)

    @abstractmethod
    def transfer_field(
        self, adata_source: AnnData, adata_target: AnnData, **kwargs
    ) -> None:
        pass

    def data_registry_mapping(self) -> dict:
        return {
            self.scvi_key: {
                _constants._DR_ATTR_NAME: self.attr_name,
                _constants._DR_ATTR_KEY: self.attr_key,
            }
        }

    def get_field(self, adata: AnnData) -> np.ndarray:
        return _get_field(adata, self.attr_name, self.attr_key)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        return dict()
