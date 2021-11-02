from __future__ import annotations

from copy import deepcopy
from typing import List, Optional, Type
from uuid import UUID, uuid4

from anndata import AnnData

import scvi

from . import _constants
from ._fields import BaseAnnDataField
from ._utils import _register_anndata, _verify_and_correct_data_format


class AnnDataManager:
    def __init__(self, fields: Optional[List[Type[BaseAnnDataField]]] = None) -> None:
        self.fields = fields or []
        self.adata = None
        self.setup_dict_key = _constants._SETUP_DICT_KEY

    def _init_setup_dict(self) -> dict:
        self._assert_anndata_registered()

        self.adata.uns[self.setup_dict_key] = {"scvi_version": scvi.__version__}

    def _assert_anndata_registered(self):
        assert (
            self.adata is not None
        ), "AnnData object not registered. Please call register_fields."

    @staticmethod
    def _validate_anndata_object(adata: AnnData):
        if adata.is_view:
            raise ValueError("Please run `adata = adata.copy()`")

    def _assign_uuid(self):
        self._assert_anndata_registered()

        if not hasattr(self.adata.uns, _constants._SCVI_UUID_KEY):
            self.adata.uns[_constants._SCVI_UUID_KEY] = uuid4()

    def _freeze_fields(self):
        self.fields = tuple(self.fields)

    def add_field(self, field: Type[BaseAnnDataField]) -> None:
        assert isinstance(
            self.fields, list
        ), "Fields have been frozen. Create a new AnnDataManager object for additional fields."
        self.fields.append(field)

    def _register_fields(
        self,
        adata: AnnData,
        source_setup_dict: Optional[dict] = None,
        **transfer_kwargs
    ):
        assert (
            self.adata is None
        ), "Existing AnnData object registered with this Manager instance."

        self._validate_anndata_object(adata)
        self.adata = adata

        self._init_setup_dict()

        for field in self.fields:
            if source_setup_dict is not None:
                field.transfer_field(source_setup_dict, self.adata, **transfer_kwargs)
            else:
                field.register_field(self.adata)
        self._freeze_fields()

        data_registry = self.get_data_registry(update=True)
        _verify_and_correct_data_format(self.adata, data_registry)

        self.get_summary_stats(update=True)

        self._assign_uuid()

    def register_fields(self, adata: AnnData):
        return self._register_fields(adata)

    def transfer_setup(
        self, adata_target: AnnData, setup_dict: Optional[dict] = None, **kwargs
    ) -> AnnDataManager:
        assert setup_dict is not None or self._assert_anndata_registered()

        setup_dict = setup_dict if setup_dict is not None else self.get_setup_dict()
        fields = self.fields
        new_adata_manager = self.__class__(fields)
        new_adata_manager._register_fields(adata_target, setup_dict, **kwargs)
        return new_adata_manager

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

        summary_stats_dict = dict()
        for field in self.fields:
            summary_stats_dict.update(field.compute_summary_stats(self.adata))

        if update:
            self.get_setup_dict()[_constants._SUMMARY_STATS_KEY] = summary_stats_dict

        return deepcopy(summary_stats_dict)
