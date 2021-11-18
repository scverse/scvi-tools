from __future__ import annotations

from copy import deepcopy
from typing import Optional, Sequence, Type
from uuid import UUID, uuid4

from anndata import AnnData

import scvi

from . import _constants
from ._utils import _register_anndata, _verify_and_correct_data_format
from .fields import BaseAnnDataField


class AnnDataManager:
    """
    Provides an interface to validate and process an AnnData object for use in scvi-tools.

    A class which wraps a collection of AnnDataField instances and provides an interface
    to validate and process an AnnData object with respect to the fields.

    Parameters
    ----------
    fields
        List of AnnDataFields to intialize with. Additional fields can be added
        via the method `add_field`.
    """

    def __init__(
        self, fields: Optional[Sequence[Type[BaseAnnDataField]]] = None
    ) -> None:
        self.fields = set(fields or {})
        self.adata = None
        self.setup_dict_key = _constants._SETUP_DICT_KEY

    def _init_setup_dict(self) -> dict:
        """Creates a setup dictionary and intializes it with the current scvi-tools version."""
        self._assert_anndata_registered()

        self.adata.uns[self.setup_dict_key] = {"scvi_version": scvi.__version__}

    def _assert_anndata_registered(self):
        """Asserts that an AnnData object has been registered with this instance."""
        assert (
            self.adata is not None
        ), "AnnData object not registered. Please call register_fields."

    @staticmethod
    def _validate_anndata_object(adata: AnnData):
        """For a given AnnData object, runs general scvi-tools compatibility checks."""
        if adata.is_view:
            raise ValueError("Please run `adata = adata.copy()`")

    def _assign_uuid(self):
        """Assigns a UUID unique to the AnnData object. If already present, the UUID is left alone."""
        self._assert_anndata_registered()

        if _constants._SCVI_UUID_KEY not in self.adata.uns:
            self.adata.uns[_constants._SCVI_UUID_KEY] = uuid4()

    def _freeze_fields(self):
        """Freezes the fields associated with this instance."""
        self.fields = frozenset(self.fields)

    def add_field(self, field: Type[BaseAnnDataField]) -> None:
        """Adds a field to this instance."""
        assert isinstance(
            self.fields, set
        ), "Fields have been frozen. Create a new AnnDataManager object for additional fields."
        self.fields.add(field)

    def _register_fields(
        self,
        adata: AnnData,
        source_setup_dict: Optional[dict] = None,
        **transfer_kwargs
    ):
        """
        Helper function with registers each field associated with this instance.

        Either registers or transfers the setup from `source_setup_dict` if passed in.

        Parameters
        ----------
        adata
            AnnData object to be registered.
        source_setup_dict
            Setup dictionary created after registering an AnnData using an AnnDataManager object.
        """
        assert (
            self.adata is None
        ), "Existing AnnData object registered with this Manager instance."

        self._validate_anndata_object(adata)
        self.adata = adata

        self._init_setup_dict()

        for field in self.fields:
            if not field.is_empty:
                if source_setup_dict is not None:
                    field.transfer_field(
                        source_setup_dict, self.adata, **transfer_kwargs
                    )
                else:
                    field.register_field(self.adata)
        self._freeze_fields()

        data_registry = self.get_data_registry(update=True)
        _verify_and_correct_data_format(self.adata, data_registry)

        self.get_summary_stats(update=True)

        self._assign_uuid()

    def register_fields(self, adata: AnnData):
        """Registers each field associated with this instance with the AnnData object."""
        return self._register_fields(adata)

    def transfer_setup(
        self, adata_target: AnnData, source_setup_dict: Optional[dict] = None, **kwargs
    ) -> AnnDataManager:
        """
        Transfers an existing setup to each field associated with this instance with the target AnnData object.

        Transfers the setup from `source_setup_dict` if passed in, otherwise uses the setup dictionary
        from the AnnData registered with this instance.

        Parameters
        ----------
        adata
            AnnData object to be registered.
        source_setup_dict
            Setup dictionary created after registering an AnnData using an AnnDataManager object.
        """
        assert source_setup_dict is not None or self.adata is not None

        setup_dict = (
            source_setup_dict
            if source_setup_dict is not None
            else self.get_setup_dict()
        )
        fields = self.fields
        new_adata_manager = self.__class__(fields)
        new_adata_manager._register_fields(adata_target, setup_dict, **kwargs)
        return new_adata_manager

    def get_adata_uuid(self) -> UUID:
        """Returns the UUID for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        return self.adata.uns[_constants._SCVI_UUID_KEY]

    def get_setup_dict(self) -> dict:
        """Returns the setup dictionary for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        return self.adata.uns[self.setup_dict_key]

    def get_data_registry(self, update: bool = False) -> dict:
        """
        Returns the data registry for the AnnData object registered with this instance.

        If `update == True`, compiles the data registry mapping and updates the data registry
        on the AnnData object.

        Parameters
        ----------
        update
            If True, recomputes and updates the data registry on the AnnData object
            registered with this instance.
        """
        self._assert_anndata_registered()

        if not update and _constants._SETUP_DICT_KEY in self.adata.uns:
            return self.get_setup_dict()

        data_registry_dict = dict()
        for field in self.fields:
            if not field.is_empty:
                data_registry_dict.update(field.data_registry_mapping())

        if update:
            _register_anndata(self.adata, data_registry_dict)

        return deepcopy(data_registry_dict)

    def get_summary_stats(self, update: bool = False) -> dict:
        """
        Returns the summary stats for the AnnData object registered with this instance.

        If `update == True`, compiles the summary stats and updates the AnnData object.

        Parameters
        ----------
        update
            If True, recomputes and updates the summary stats on the AnnData object
            registered with this instance.
        """
        self._assert_anndata_registered()

        setup_dict = self.get_setup_dict()

        if not update and _constants._SUMMARY_STATS_KEY in setup_dict:
            return self.get_setup_dict()[_constants._SUMMARY_STATS_KEY]

        summary_stats_dict = dict()
        for field in self.fields:
            summary_stats_dict.update(field.compute_summary_stats(self.adata))

        if update:
            setup_dict[_constants._SUMMARY_STATS_KEY] = summary_stats_dict

        return deepcopy(summary_stats_dict)
