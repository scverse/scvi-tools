from copy import deepcopy
from typing import Type

from anndata import AnnData

import scvi

from . import _constants
from ._fields import BaseAnnDataField
from ._utils import _register_anndata, _verify_and_correct_data_format


class AnnDataManager:
    def __init__(self, adata: AnnData) -> None:
        self.fields = []
        self.adata = adata

        self._init_setup_dict()

    def _init_setup_dict(self) -> dict:
        self.adata.uns[_constants._SETUP_DICT_KEY] = {"scvi_version": scvi.__version__}

    def add_field(self, field: Type[BaseAnnDataField]) -> None:
        self.fields.append(field)

    def get_data_registry(self, update: bool = False) -> dict:
        data_registry_dict = dict()
        for field in self.fields:
            data_registry_dict.update(field.data_registry_mapping())

        if update:
            _register_anndata(self.adata, data_registry_dict)

        return deepcopy(data_registry_dict)

    def get_summary_stats(self, update: bool = False) -> dict:
        n_cells, n_vars = self.adata.shape
        summary_stats_dict = dict(n_cells=n_cells, n_vars=n_vars)
        for field in self.fields:
            summary_stats_dict.update(field.compute_summary_stats(self.adata))

        if update:
            self.adata.uns[_constants._SETUP_DICT_KEY][
                _constants._SUMMARY_STATS_KEY
            ] = summary_stats_dict

        return deepcopy(summary_stats_dict)

    def register_fields(self) -> None:
        for field in self.fields:
            field.register_field(self.adata)

        data_registry = self.get_data_registry(update=True)
        _verify_and_correct_data_format(self.adata, data_registry)

        self.get_summary_stats(update=True)
