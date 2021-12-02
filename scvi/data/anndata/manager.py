from __future__ import annotations

from copy import deepcopy
from typing import Optional, Sequence, Type
from uuid import UUID, uuid4

from anndata import AnnData

import scvi

from . import _constants
from ._utils import _verify_and_correct_data_format
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
        self.registry = {_constants._SCVI_VERSION_KEY: scvi.__version__}

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

        scvi_uuid = self.adata.uns[_constants._SCVI_UUID_KEY]
        self.registry[_constants._SCVI_UUID_KEY] = scvi_uuid

    def _freeze_fields(self):
        """Freezes the fields associated with this instance."""
        self.fields = frozenset(self.fields)

    def add_field(self, field: Type[BaseAnnDataField]) -> None:
        """Adds a field to this instance."""
        assert isinstance(
            self.fields, set
        ), "Fields have been frozen. Create a new AnnDataManager object for additional fields."
        self.fields.add(field)

    def register_fields(
        self, adata: AnnData, source_registry: Optional[dict] = None, **transfer_kwargs
    ):
        """
        Registers each field associated with this instance with the AnnData object.

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

        for field in self.fields:
            self.registry[field.registry_key] = {
                _constants._DATA_REGISTRY_KEY: field.get_data_registry(),
                _constants._STATE_REGISTRY_KEY: dict(),
            }
            field_registry = self.registry[field.registry_key]

            if not field.is_empty:
                if source_registry is not None:
                    field_registry[
                        _constants._STATE_REGISTRY_KEY
                    ] = field.transfer_field(
                        source_registry[field.registry_key][
                            _constants._STATE_REGISTRY_KEY
                        ],
                        self.adata,
                        **transfer_kwargs
                    )
                else:
                    field_registry[
                        _constants._STATE_REGISTRY_KEY
                    ] = field.register_field(self.adata)

            state_registry = field_registry[_constants._STATE_REGISTRY_KEY]
            field_registry[_constants._SUMMARY_STATS_KEY] = field.get_summary_stats(
                state_registry
            )

        self._freeze_fields()

        _verify_and_correct_data_format(self.adata, self.data_registry)

        self._assign_uuid()

    def transfer_setup(
        self, adata_target: AnnData, source_registry: Optional[dict] = None, **kwargs
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
        assert source_registry is not None or self.adata is not None

        fields = self.fields
        new_adata_manager = self.__class__(fields)
        new_adata_manager.register_fields(adata_target, source_registry, **kwargs)
        return new_adata_manager

    def get_adata_uuid(self) -> UUID:
        """Returns the UUID for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        return self.adata.uns[_constants._SCVI_UUID_KEY]

    @property
    def data_registry(self) -> dict:
        """Returns the data registry for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        data_registry = dict()
        for registry_key, field_registry in self.registry:
            field_data_registry = field_registry[_constants._DATA_REGISTRY_KEY]
            if field_data_registry:
                data_registry[registry_key] = field_data_registry.copy()

        return data_registry

    @property
    def summary_stats(self) -> dict:
        """Returns the summary stats for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        summary_stats = dict()
        for field_registry in self.registry.values():
            field_summary_stats = field_registry[_constants._SUMMARY_STATS_KEY]
            summary_stats.update(field_summary_stats)

        return deepcopy(summary_stats)
