from __future__ import annotations

import sys
from collections import defaultdict
from copy import deepcopy
from typing import Optional, Sequence, Type, Union
from uuid import uuid4

import numpy as np
import pandas as pd
import rich
from anndata import AnnData

import scvi
from scvi.utils import attrdict

from . import _constants
from ._utils import _assign_adata_uuid, get_anndata_attribute
from .fields import BaseAnnDataField


class AnnDataManager:
    """
    Provides an interface to validate and process an AnnData object for use in scvi-tools.

    A class which wraps a collection of AnnDataField instances and provides an interface
    to validate and process an AnnData object with respect to the fields.

    Parameters
    ----------
    fields
        List of AnnDataFields to intialize with.
    setup_method_args
        Dictionary describing the model and arguments passed in by the user
        to setup this AnnDataManager.

    Examples
    --------
    >>> fields = [LayerField("counts", "raw_counts")]
    >>> adata_manager = AnnDataManager(fields=fields)
    >>> adata_manager.register_fields(adata)

    Notes
    -----
    This class is not initialized with a specific AnnData object, but later sets ``self.adata``
    via :meth:`~scvi.data.AnnDataManager.register_fields`. This decouples the generalized
    definition of the scvi-tools interface with the registration of an instance of data.
    """

    def __init__(
        self,
        fields: Optional[Sequence[Type[BaseAnnDataField]]] = None,
        setup_method_args: Optional[dict] = None,
    ) -> None:
        self.id = str(uuid4())
        self.adata = None
        self.fields = fields or []
        self._registry = {
            _constants._SCVI_VERSION_KEY: scvi.__version__,
            _constants._MODEL_NAME_KEY: None,
            _constants._SETUP_ARGS_KEY: None,
            _constants._FIELD_REGISTRIES_KEY: defaultdict(dict),
        }
        if setup_method_args is not None:
            self._registry.update(setup_method_args)

    def _assert_anndata_registered(self):
        """Asserts that an AnnData object has been registered with this instance."""
        if self.adata is None:
            raise AssertionError(
                "AnnData object not registered. Please call register_fields."
            )

    @staticmethod
    def _validate_anndata_object(adata: AnnData):
        """For a given AnnData object, runs general scvi-tools compatibility checks."""
        if adata.is_view:
            raise ValueError("Please run `adata = adata.copy()`")

    def _get_setup_method_args(self) -> dict:
        """
        Returns the ``setup_anndata`` method arguments used to initialize this :class:`~scvi.data.AnnDataManager` instance.

        Returns the ``setup_anndata`` method arguments, including the model name,
        that were used to initialize this :class:`~scvi.data.AnnDataManager` instance
        in the form of a dictionary.
        """
        return {
            k: v
            for k, v in self._registry.items()
            if k in {_constants._MODEL_NAME_KEY, _constants._SETUP_ARGS_KEY}
        }

    def _assign_uuid(self):
        """Assigns a UUID unique to the AnnData object. If already present, the UUID is left alone."""
        self._assert_anndata_registered()

        _assign_adata_uuid(self.adata)

        scvi_uuid = self.adata.uns[_constants._SCVI_UUID_KEY]
        self._registry[_constants._SCVI_UUID_KEY] = scvi_uuid

    def _assign_most_recent_manager_uuid(self):
        """
        Assigns a last manager UUID to the AnnData object for future validation.
        """
        self._assert_anndata_registered()

        self.adata.uns[_constants._MANAGER_UUID_KEY] = self.id

    def register_fields(
        self, adata: AnnData, source_registry: Optional[dict] = None, **transfer_kwargs
    ):
        """
        Registers each field associated with this instance with the AnnData object.

        Either registers or transfers the setup from `source_setup_dict` if passed in.
        Sets ``self.adata``.

        Parameters
        ----------
        adata
            AnnData object to be registered.
        source_registry
            Registry created after registering an AnnData using an :class:`~scvi.data.AnnDataManager` object.
        transfer_kwargs
            Additional keywords which modify transfer behavior. Only applicable if ``source_registry`` is set.
        """
        if self.adata is not None:
            raise AssertionError(
                "Existing AnnData object registered with this Manager instance."
            )

        if source_registry is None and transfer_kwargs:
            raise TypeError(
                f"register_fields() got unexpected keyword arguments {transfer_kwargs} passed without a source_registry."
            )

        self._validate_anndata_object(adata)
        field_registries = self._registry[_constants._FIELD_REGISTRIES_KEY]

        for field in self.fields:
            field_registries[field.registry_key] = {
                _constants._DATA_REGISTRY_KEY: field.get_data_registry(),
                _constants._STATE_REGISTRY_KEY: dict(),
            }
            field_registry = field_registries[field.registry_key]

            # A field can be empty if the model has optional fields (e.g. extra covariates).
            # If empty, we skip registering the field.
            if not field.is_empty:
                # Transfer case: Source registry is used for validation and/or setup.
                if source_registry is not None:
                    field_registry[
                        _constants._STATE_REGISTRY_KEY
                    ] = field.transfer_field(
                        source_registry[_constants._FIELD_REGISTRIES_KEY][
                            field.registry_key
                        ][_constants._STATE_REGISTRY_KEY],
                        adata,
                        **transfer_kwargs,
                    )
                else:
                    field_registry[
                        _constants._STATE_REGISTRY_KEY
                    ] = field.register_field(adata)

            # Compute and set summary stats for the given field.
            state_registry = field_registry[_constants._STATE_REGISTRY_KEY]
            field_registry[_constants._SUMMARY_STATS_KEY] = field.get_summary_stats(
                state_registry
            )

        # Save arguments for register_fields.
        self._source_registry = deepcopy(source_registry)
        self._transfer_kwargs = deepcopy(transfer_kwargs)

        self.adata = adata
        self._assign_uuid()
        self._assign_most_recent_manager_uuid()

    def transfer_fields(self, adata_target: AnnData, **kwargs) -> AnnDataManager:
        """
        Transfers an existing setup to each field associated with this instance with the target AnnData object.

        Creates a new :class:`~scvi.data.AnnDataManager` instance with the same set of fields.
        Then, registers the fields with a target AnnData object, incorporating details of the
        source registry where necessary (e.g. for validation or modified data setup).

        Parameters
        ----------
        adata_target
            AnnData object to be registered.
        kwargs
            Additional keywords which modify transfer behavior.
        """
        self._assert_anndata_registered()

        fields = self.fields
        new_adata_manager = self.__class__(
            fields=fields, setup_method_args=self._get_setup_method_args()
        )
        new_adata_manager.register_fields(adata_target, self._registry, **kwargs)
        return new_adata_manager

    def validate(self) -> None:
        """Checks if AnnData was last setup with this AnnDataManager instance and reregisters it if not."""
        self._assert_anndata_registered()
        most_recent_manager_id = self.adata.uns[_constants._MANAGER_UUID_KEY]
        # Re-register fields with same arguments if this AnnData object has been
        # registered with a different AnnDataManager.
        if most_recent_manager_id != self.id:
            adata, self.adata = self.adata, None  # Reset self.adata.
            self.register_fields(adata, self._source_registry, **self._transfer_kwargs)

    @property
    def adata_uuid(self) -> str:
        """Returns the UUID for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        return self._registry[_constants._SCVI_UUID_KEY]

    @property
    def registry(self) -> dict:
        """Returns the top-level registry dictionary for the AnnData object registered with this instance as an attrdict."""
        return self._registry

    @property
    def data_registry(self) -> attrdict:
        """Returns the data registry for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        data_registry = dict()
        for registry_key, field_registry in self._registry[
            _constants._FIELD_REGISTRIES_KEY
        ].items():
            field_data_registry = field_registry[_constants._DATA_REGISTRY_KEY]
            if field_data_registry:
                data_registry[registry_key] = field_data_registry

        return attrdict(data_registry, recursive=True)

    @property
    def summary_stats(self) -> attrdict:
        """Returns the summary stats for the AnnData object registered with this instance."""
        self._assert_anndata_registered()

        summary_stats = dict()
        for field_registry in self._registry[_constants._FIELD_REGISTRIES_KEY].values():
            field_summary_stats = field_registry[_constants._SUMMARY_STATS_KEY]
            summary_stats.update(field_summary_stats)

        return attrdict(summary_stats)

    def get_from_registry(self, registry_key: str) -> Union[np.ndarray, pd.DataFrame]:
        """
        Returns the object in AnnData associated with the key in the data registry.

        Parameters
        ----------
        registry_key
            key of object to get from ``self.data_registry``

        Returns
        -------
        The requested data.
        """
        data_loc = self.data_registry[registry_key]
        attr_name, attr_key = (
            data_loc[_constants._DR_ATTR_NAME],
            data_loc[_constants._DR_ATTR_KEY],
        )

        return get_anndata_attribute(self.adata, attr_name, attr_key)

    def get_state_registry(self, registry_key: str) -> attrdict:
        """Returns the state registry for the AnnDataField registered with this instance."""
        self._assert_anndata_registered()

        return attrdict(
            self._registry[_constants._FIELD_REGISTRIES_KEY][registry_key][
                _constants._STATE_REGISTRY_KEY
            ]
        )

    def _view_summary_stats(self) -> rich.table.Table:
        """Prints summary stats."""
        t = rich.table.Table(title="Summary Statistics")
        t.add_column(
            "Summary Stat Key",
            justify="center",
            style="dodger_blue1",
            no_wrap=True,
            overflow="fold",
        )
        t.add_column(
            "Value",
            justify="center",
            style="dark_violet",
            no_wrap=True,
            overflow="fold",
        )
        for stat_key, count in self.summary_stats.items():
            t.add_row(stat_key, str(count))
        return t

    def _view_data_registry(self) -> rich.table.Table:
        """Prints data registry."""
        t = rich.table.Table(title="Data Registry")
        t.add_column(
            "Registry Key",
            justify="center",
            style="dodger_blue1",
            no_wrap=True,
            overflow="fold",
        )
        t.add_column(
            "scvi-tools Location",
            justify="center",
            style="dark_violet",
            no_wrap=True,
            overflow="fold",
        )

        for registry_key, data_loc in self.data_registry.items():
            attr_name = data_loc.attr_name
            attr_key = data_loc.attr_key
            if attr_key is None:
                scvi_data_str = f"adata.{attr_name}"
            else:
                scvi_data_str = f"adata.{attr_name}['{attr_key}']"
            t.add_row(registry_key, scvi_data_str)

        return t

    @staticmethod
    def view_setup_method_args(registry: dict) -> None:
        """
        Prints setup kwargs used to produce a given registry.

        Parameters
        ----------
        registry
            Registry produced by an AnnDataManager.
        """
        model_name = registry[_constants._MODEL_NAME_KEY]
        setup_args = registry[_constants._SETUP_ARGS_KEY]
        if model_name is not None and setup_args is not None:
            rich.print(f"Setup via `{model_name}.setup_anndata` with arguments:")
            rich.pretty.pprint(setup_args)
            rich.print()

    def view_registry(self, hide_state_registries: bool = False) -> None:
        """
        Prints summary of the registry.

        Parameters
        ----------
        hide_state_registries
            If True, prints a shortened summary without details of each state registry.
        """
        version = self._registry[_constants._SCVI_VERSION_KEY]
        rich.print(f"Anndata setup with scvi-tools version {version}.")
        rich.print()
        self.view_setup_method_args(self._registry)

        in_colab = "google.colab" in sys.modules
        force_jupyter = None if not in_colab else True
        console = rich.console.Console(force_jupyter=force_jupyter)
        console.print(self._view_summary_stats())
        console.print(self._view_data_registry())

        if not hide_state_registries:
            for field in self.fields:
                state_registry = self.get_state_registry(field.registry_key)
                t = field.view_state_registry(state_registry)
                if t is not None:
                    console.print(t)
