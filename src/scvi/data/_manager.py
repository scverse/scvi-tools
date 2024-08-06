from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass
from uuid import uuid4

import numpy as np
import pandas as pd
from mudata import MuData
from torch.utils.data import Subset

import scvi
from scvi._types import AnnOrMuData
from scvi.utils import attrdict

from . import _constants
from ._anntorchdataset import AnnTorchDataset
from ._utils import (
    _assign_adata_uuid,
    _check_if_view,
    _check_mudata_fully_paired,
    get_anndata_attribute,
)
from .fields import AnnDataField


@dataclass
class AnnDataManagerValidationCheck:
    """Validation checks for AnnorMudata scvi-tools compat.

    Parameters
    ----------
    check_if_view
        If True, checks if AnnData is a view.
    check_fully_paired_mudata
        If True, checks if MuData is fully paired across mods.
    """

    check_if_view: bool = True
    check_fully_paired_mudata: bool = True


class AnnDataManager:
    """Provides an interface to validate and process an AnnData object for use in scvi-tools.

    A class which wraps a collection of AnnDataField instances and provides an interface
    to validate and process an AnnData object with respect to the fields.

    Parameters
    ----------
    fields
        List of AnnDataFields to intialize with.
    setup_method_args
        Dictionary describing the model and arguments passed in by the user
        to setup this AnnDataManager.
    validation_checks
        DataClass specifying which global validation checks to run on the data object.

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

    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/dev/data_tutorial`
    """

    def __init__(
        self,
        fields: list[AnnDataField] | None = None,
        setup_method_args: dict | None = None,
        validation_checks: AnnDataManagerValidationCheck | None = None,
    ) -> None:
        self.id = str(uuid4())
        self.adata = None
        self.fields = fields or []
        self.validation_checks = validation_checks or AnnDataManagerValidationCheck()
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
            raise AssertionError("AnnData object not registered. Please call register_fields.")

    def _validate_anndata_object(self, adata: AnnOrMuData):
        """For a given AnnData object, runs general scvi-tools compatibility checks."""
        if self.validation_checks.check_if_view:
            _check_if_view(adata, copy_if_view=False)

        if isinstance(adata, MuData) and self.validation_checks.check_fully_paired_mudata:
            _check_mudata_fully_paired(adata)

    def _get_setup_method_args(self) -> dict:
        """Returns the ``setup_anndata`` method arguments.

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
        """Assigns a UUID unique to the AnnData object.

        If already present, the UUID is left alone.
        """
        self._assert_anndata_registered()

        _assign_adata_uuid(self.adata)

        scvi_uuid = self.adata.uns[_constants._SCVI_UUID_KEY]
        self._registry[_constants._SCVI_UUID_KEY] = scvi_uuid

    def _assign_most_recent_manager_uuid(self):
        """Assigns a last manager UUID to the AnnData object for future validation."""
        self._assert_anndata_registered()

        self.adata.uns[_constants._MANAGER_UUID_KEY] = self.id

    def register_fields(
        self,
        adata: AnnOrMuData,
        source_registry: dict | None = None,
        **transfer_kwargs,
    ):
        """Registers each field associated with this instance with the AnnData object.

        Either registers or transfers the setup from `source_setup_dict` if passed in.
        Sets ``self.adata``.

        Parameters
        ----------
        adata
            AnnData object to be registered.
        source_registry
            Registry created after registering an AnnData using an
            :class:`~scvi.data.AnnDataManager` object.
        transfer_kwargs
            Additional keywords which modify transfer behavior. Only applicable if
            ``source_registry`` is set.
        """
        if self.adata is not None:
            raise AssertionError("Existing AnnData object registered with this Manager instance.")

        if source_registry is None and transfer_kwargs:
            raise TypeError(
                f"register_fields() got unexpected keyword arguments {transfer_kwargs} passed "
                "without a source_registry."
            )

        self._validate_anndata_object(adata)

        for field in self.fields:
            self._add_field(
                field=field,
                adata=adata,
                source_registry=source_registry,
                **transfer_kwargs,
            )

        # Save arguments for register_fields.
        self._source_registry = deepcopy(source_registry)
        self._transfer_kwargs = deepcopy(transfer_kwargs)

        self.adata = adata
        self._assign_uuid()
        self._assign_most_recent_manager_uuid()

    def _add_field(
        self,
        field: AnnDataField,
        adata: AnnOrMuData,
        source_registry: dict | None = None,
        **transfer_kwargs,
    ):
        """Internal function for adding a field with optional transferring."""
        field_registries = self._registry[_constants._FIELD_REGISTRIES_KEY]
        field_registries[field.registry_key] = {
            _constants._DATA_REGISTRY_KEY: field.get_data_registry(),
            _constants._STATE_REGISTRY_KEY: {},
        }
        field_registry = field_registries[field.registry_key]

        # A field can be empty if the model has optional fields (e.g. extra covariates).
        # If empty, we skip registering the field.
        if not field.is_empty:
            # Transfer case: Source registry is used for validation and/or setup.
            if source_registry is not None:
                field_registry[_constants._STATE_REGISTRY_KEY] = field.transfer_field(
                    source_registry[_constants._FIELD_REGISTRIES_KEY][field.registry_key][
                        _constants._STATE_REGISTRY_KEY
                    ],
                    adata,
                    **transfer_kwargs,
                )
            else:
                field_registry[_constants._STATE_REGISTRY_KEY] = field.register_field(adata)
        # Compute and set summary stats for the given field.
        state_registry = field_registry[_constants._STATE_REGISTRY_KEY]
        field_registry[_constants._SUMMARY_STATS_KEY] = field.get_summary_stats(state_registry)

    def register_new_fields(self, fields: list[AnnDataField]):
        """Register new fields to a manager instance.

        This is useful to augment the functionality of an existing manager.

        Parameters
        ----------
        fields
            List of AnnDataFields to register
        """
        if self.adata is None:
            raise AssertionError(
                "No AnnData object has been registered with this Manager instance."
            )
        self.validate()
        for field in fields:
            self._add_field(
                field=field,
                adata=self.adata,
            )

        # Source registry is not None if this manager was created from transfer_fields
        # In this case self._registry is originally equivalent to self._source_registry
        # However, with newly registered fields the equality breaks so we reset it
        if self._source_registry is not None:
            self._source_registry = deepcopy(self._registry)

        self.fields += fields

    def transfer_fields(self, adata_target: AnnOrMuData, **kwargs) -> AnnDataManager:
        """Transfers an existing setup to each field associated with this instance.

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
            fields=fields,
            setup_method_args=self._get_setup_method_args(),
            validation_checks=self.validation_checks,
        )
        new_adata_manager.register_fields(adata_target, self._registry, **kwargs)
        return new_adata_manager

    def validate(self) -> None:
        """Checks if AnnData was last setup with this AnnDataManager instanc.

        Reregisters it if not.
        """
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
    def data_registry(self) -> attrdict:
        """Returns the data registry for the AnnData object registered with this instance."""
        self._assert_anndata_registered()
        return self._get_data_registry_from_registry(self._registry)

    def create_torch_dataset(
        self,
        indices: Sequence[int] | Sequence[bool] = None,
        data_and_attributes: list[str] | dict[str, np.dtype] | None = None,
        load_sparse_tensor: bool = False,
    ) -> AnnTorchDataset:
        """
        Creates a torch dataset from the AnnData object registered with this instance.

        Parameters
        ----------
        indices
            The indices of the observations in the adata to use
        data_and_attributes
            Dictionary with keys representing keys in data registry
            (``adata_manager.data_registry``) and value equal to desired numpy loading type (later
            made into torch tensor) or list of such keys. A list can be used to subset to certain
            keys in the event that more tensors than needed have been registered. If ``None``,
            defaults to all registered data.
        load_sparse_tensor
            ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
            :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to
            GPUs, depending on the sparsity of the data.

        Returns
        -------
        :class:`~scvi.data.AnnTorchDataset`
        """
        dataset = AnnTorchDataset(
            self,
            getitem_tensors=data_and_attributes,
            load_sparse_tensor=load_sparse_tensor,
        )
        if indices is not None:
            # This is a lazy subset, it just remaps indices
            dataset = Subset(dataset, indices)
        return dataset

    @staticmethod
    def _get_data_registry_from_registry(registry: dict) -> attrdict:
        data_registry = {}
        for registry_key, field_registry in registry[_constants._FIELD_REGISTRIES_KEY].items():
            field_data_registry = field_registry[_constants._DATA_REGISTRY_KEY]
            if field_data_registry:
                data_registry[registry_key] = field_data_registry
        return attrdict(data_registry)

    def get_from_registry(self, registry_key: str) -> np.ndarray | pd.DataFrame:
        """Returns the object in AnnData associated with the key in the data registry.

        Parameters
        ----------
        registry_key
            key of object to get from ``self.data_registry``

        Returns
        -------
        The requested data.
        """
        data_loc = self.data_registry[registry_key]
        mod_key, attr_name, attr_key = (
            getattr(data_loc, _constants._DR_MOD_KEY, None),
            data_loc[_constants._DR_ATTR_NAME],
            data_loc[_constants._DR_ATTR_KEY],
        )

        return get_anndata_attribute(self.adata, attr_name, attr_key, mod_key=mod_key)

