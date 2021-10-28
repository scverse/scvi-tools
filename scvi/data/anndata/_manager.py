from copy import deepcopy
from typing import Type
from uuid import UUID, uuid4

from anndata import AnnData

import scvi

from . import _constants
from ._fields import BaseAnnDataField
from ._utils import _register_anndata, _verify_and_correct_data_format


class AnnDataManager:
    def __init__(self) -> None:
        self.fields = []
        self.adata = None
        self.setup_dict_key = _constants._SETUP_DICT_KEY

    def _init_setup_dict(self) -> dict:
        self._assert_anndata_registered()

        self.adata.uns[self.setup_dict_key] = {"scvi_version": scvi.__version__}

    def _assert_anndata_registered(self):
        assert (
            self.adata is not None
        ), "AnnData object not registered. Please call register_fields."

    def _add_field(self, field: Type[BaseAnnDataField]) -> None:
        self.fields.append(field)

    @staticmethod
    def _validate_anndata_object(adata: AnnData):
        if adata.is_view:
            raise ValueError("Please run `adata = adata.copy()`")

    def _assign_uuid(self):
        self._assert_anndata_registered()

        if not hasattr(self.adata.uns, _constants._SCVI_UUID_KEY):
            self.adata.uns[_constants._SCVI_UUID_KEY] = uuid4()

    def _register_fields(self, adata: AnnData):
        assert (
            self.adata is None
        ), "Existing AnnData object registered with this Manager instance."

        self._validate_anndata_object(adata)
        self.adata = adata

        self._init_setup_dict()

        for field in self.fields:
            field.register_field(self.adata)

        data_registry = self.get_data_registry(update=True)
        _verify_and_correct_data_format(self.adata, data_registry)

        self.get_summary_stats(update=True)

        self._assign_uuid()

    def get_adata_uuid(self) -> UUID:
        self._assert_anndata_registered()

        return self.adata.uns[_constants._SCVI_UUID_KEY]

    def get_setup_dict(self) -> dict:
        self._assert_anndata_registered()

        return self.adata.uns[self.setup_dict_key]

    def get_data_registry(self, update: bool = False) -> dict:
        self._assert_anndata_registered()

        data_registry_dict = dict()
        for field in self.fields:
            data_registry_dict.update(field.data_registry_mapping())

        if update:
            _register_anndata(self.adata, data_registry_dict)

        return deepcopy(data_registry_dict)

    def get_summary_stats(self, update: bool = False) -> dict:
        self._assert_anndata_registered()

        n_cells, n_vars = self.adata.shape
        summary_stats_dict = dict(n_cells=n_cells, n_vars=n_vars)
        for field in self.fields:
            summary_stats_dict.update(field.compute_summary_stats(self.adata))

        if update:
            self.get_setup_dict()[_constants._SUMMARY_STATS_KEY] = summary_stats_dict

        return deepcopy(summary_stats_dict)
