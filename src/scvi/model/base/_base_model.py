from __future__ import annotations

import inspect
import logging
import os
import sys
import warnings
from abc import ABCMeta, abstractmethod
from io import StringIO
from typing import TYPE_CHECKING
from uuid import uuid4

import numpy as np
import pyro
import rich
import torch
from anndata import AnnData
from mudata import MuData
from rich import box
from rich.console import Console

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager, fields
from scvi.data._compat import registry_from_setup_dict
from scvi.data._constants import (
    _ADATA_MINIFY_TYPE_UNS_KEY,
    _FIELD_REGISTRIES_KEY,
    _MODEL_NAME_KEY,
    _SCVI_UUID_KEY,
    _SCVI_VERSION_KEY,
    _SETUP_ARGS_KEY,
    _SETUP_METHOD_NAME,
    _STATE_REGISTRY_KEY,
    ADATA_MINIFY_TYPE,
)
from scvi.data._utils import (
    _assign_adata_uuid,
    _check_if_view,
    _get_adata_minify_type,
)
from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_device_args
from scvi.model.base._constants import SAVE_KEYS
from scvi.model.base._save_load import (
    _initialize_model,
    _load_legacy_saved_files,
    _load_saved_files,
    _validate_var_names,
)
from scvi.model.utils import get_minified_adata_scrna, get_minified_mudata
from scvi.utils import attrdict, setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

from . import _constants

if TYPE_CHECKING:
    from collections.abc import Sequence

    import pandas as pd
    from lightning import LightningDataModule

    from scvi._types import AnnOrMuData, MinifiedDataType

logger = logging.getLogger(__name__)


_UNTRAINED_WARNING_MESSAGE = (
    "Trying to query inferred values from an untrained model. Please train the model first."
)

_SETUP_INPUTS_EXCLUDED_PARAMS = {"adata", "mdata", "kwargs"}


class BaseModelMetaClass(ABCMeta):
    """Metaclass for :class:`~scvi.model.base.BaseModelClass`.

    Constructs model class-specific mappings for :class:`~scvi.data.AnnDataManager` instances.
    ``cls._setup_adata_manager_store`` maps from AnnData object UUIDs to
    :class:`~scvi.data.AnnDataManager` instances.

    This mapping is populated everytime ``cls.setup_anndata()`` is called.
    ``cls._per_isntance_manager_store`` maps from model instance UUIDs to AnnData UUID:
    :class:`~scvi.data.AnnDataManager` mappings.
    These :class:`~scvi.data.AnnDataManager` instances are tied to a single model instance and
    populated either
    during model initialization or after running ``self._validate_anndata()``.
    """

    @abstractmethod
    def __init__(cls, name, bases, dct):
        cls._setup_adata_manager_store: dict[
            str, type[AnnDataManager]
        ] = {}  # Maps adata id to AnnDataManager instances.
        cls._per_instance_manager_store: dict[
            str, dict[str, type[AnnDataManager]]
        ] = {}  # Maps model instance id to AnnDataManager mappings.
        super().__init__(name, bases, dct)


class BaseModelClass(metaclass=BaseModelMetaClass):
    """Abstract class for scvi-tools models.

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/dev/model_user_guide`
    """

    _OBSERVED_LIB_SIZE_KEY = "observed_lib_size"
    _data_loader_cls = AnnDataLoader

    def __init__(self, adata: AnnOrMuData | None = None, registry: object | None = None):
        # check if the given adata is minified and check if the model being created
        # supports minified-data mode (i.e. inherits from the abstract BaseMinifiedModeModelClass).
        # If not, raise an error to inform the user of the lack of minified-data functionality
        # for this model
        data_is_minified = adata is not None and _get_adata_minify_type(adata) is not None
        if data_is_minified and not issubclass(type(self), BaseMinifiedModeModelClass):
            raise NotImplementedError(
                f"The {type(self).__name__} model currently does not support minified data."
            )
        self.id = str(uuid4())  # Used for cls._manager_store keys.
        if adata is not None:
            self._adata = adata
            self._adata_manager = self._get_most_recent_anndata_manager(adata, required=True)
            self._register_manager_for_instance(self.adata_manager)
            # Suffix registry instance variable with _ to include it when saving the model.
            self.registry_ = self._adata_manager._registry
            self.summary_stats = AnnDataManager._get_summary_stats_from_registry(self.registry_)
        elif registry is not None:
            self._adata = None
            self._adata_manager = None
            # Suffix registry instance variable with _ to include it when saving the model.
            self.registry_ = registry
            self.summary_stats = AnnDataManager._get_summary_stats_from_registry(registry)
        elif (self.__class__.__name__ == "GIMVI") or (self.__class__.__name__ == "SCVI"):
            # note some models do accept empty registry/adata (e.g: gimvi)
            pass
        else:
            raise ValueError("adata or registry must be provided.")

        self._module_init_on_train = adata is None and registry is None
        self.is_trained_ = False
        self._model_summary_string = ""
        self.train_indices_ = None
        self.test_indices_ = None
        self.validation_indices_ = None
        self.history_ = None
        self.get_normalized_function_name_ = "get_normalized_expression"

    @property
    def adata(self) -> None | AnnOrMuData:
        """Data attached to model instance."""
        return self._adata

    @property
    def registry(self) -> dict:
        """Data attached to model instance."""
        return self.registry_

    def get_var_names(self, legacy_mudata_format=False) -> dict:
        """Variable names of input data."""
        from scvi.model.base._save_load import _get_var_names

        if self.adata:
            return _get_var_names(self.adata, legacy_mudata_format=legacy_mudata_format)
        else:
            return self.registry[_FIELD_REGISTRIES_KEY]["X"][_STATE_REGISTRY_KEY]["column_names"]

    @adata.setter
    def adata(self, adata: AnnOrMuData):
        if adata is None:
            raise ValueError("adata cannot be None.")
        self._validate_anndata(adata)
        self._adata = adata
        self._adata_manager = self.get_anndata_manager(adata, required=True)
        self.registry_ = self._adata_manager.registry
        self.summary_stats = self._adata_manager.summary_stats

    @property
    def adata_manager(self) -> AnnDataManager:
        """Manager instance associated with self.adata."""
        return self._adata_manager

    def to_device(self, device: str | int):
        """Move model to device.

        Parameters
        ----------
        device
            Device to move model to. Options: 'cpu' for CPU, integer GPU index (eg. 0),
            or 'cuda:X' where X is the GPU index (eg. 'cuda:0'). See torch.device for more info.

        Examples
        --------
        >>> adata = scvi.data.synthetic_iid()
        >>> model = scvi.model.SCVI(adata)
        >>> model.to_device("cpu")  # moves model to CPU
        >>> model.to_device("cuda:0")  # moves model to GPU 0
        >>> model.to_device(0)  # also moves model to GPU 0
        """
        my_device = torch.device(device)
        self.module.to(my_device)

    @property
    def device(self) -> str:
        """The current device that the module's params are on."""
        return self.module.device

    @staticmethod
    def _get_setup_method_args(**setup_locals) -> dict:
        """Returns a dictionary organizing the arguments used to call ``setup_anndata``.

        Must be called with ``**locals()`` at the start of the ``setup_anndata`` method
        to avoid the inclusion of any extraneous variables.
        """
        cls = setup_locals.pop("cls")
        method_name = None
        if "adata" in setup_locals:
            method_name = "setup_anndata"
        elif "mdata" in setup_locals:
            method_name = "setup_mudata"

        model_name = cls.__name__
        setup_args = {}
        for k, v in setup_locals.items():
            if k not in _SETUP_INPUTS_EXCLUDED_PARAMS:
                setup_args[k] = v
        return {
            _MODEL_NAME_KEY: model_name,
            _SETUP_METHOD_NAME: method_name,
            _SETUP_ARGS_KEY: setup_args,
        }

    @staticmethod
    def _create_modalities_attr_dict(
        modalities: dict[str, str], setup_method_args: dict
    ) -> attrdict:
        """Preprocesses a ``modalities`` dictionary to map modality names.

        Ensures each field key has a respective modality key, defaulting to ``None``.
        Raises a ``UserWarning`` if extraneous modality mappings are detected.

        Parameters
        ----------
        modalities
            Dictionary mapping ``setup_mudata()`` argument name to modality name.
        setup_method_args
            Output of  ``_get_setup_method_args()``.
        """
        setup_args = setup_method_args[_SETUP_ARGS_KEY]
        filtered_modalities = {
            arg_name: modalities.get(arg_name, None) for arg_name in setup_args.keys()
        }
        extra_modalities = set(modalities) - set(filtered_modalities)
        if len(extra_modalities) > 0:
            raise ValueError(f"Extraneous modality mapping(s) detected: {extra_modalities}")
        return attrdict(filtered_modalities)

    @classmethod
    def register_manager(cls, adata_manager: AnnDataManager):
        """Registers an :class:`~scvi.data.AnnDataManager` instance with this model class.

        Stores the :class:`~scvi.data.AnnDataManager` reference in a class-specific manager store.
        Intended for use in the ``setup_anndata()`` class method followed up by retrieval of the
        :class:`~scvi.data.AnnDataManager` via the ``_get_most_recent_anndata_manager()`` method in
        the model init method.

        Notes
        -----
        Subsequent calls to this method with an :class:`~scvi.data.AnnDataManager` instance
        referring to the same underlying AnnData object will overwrite the reference to previous
        :class:`~scvi.data.AnnDataManager`.
        """
        adata_id = adata_manager.adata_uuid
        cls._setup_adata_manager_store[adata_id] = adata_manager

    def _register_manager_for_instance(self, adata_manager: AnnDataManager):
        """Registers an :class:`~scvi.data.AnnDataManager` instance with this model instance.

        Creates a model-instance specific mapping in ``cls._per_instance_manager_store`` for this
        :class:`~scvi.data.AnnDataManager` instance.
        """
        if self.id not in self._per_instance_manager_store:
            self._per_instance_manager_store[self.id] = {}

        adata_id = adata_manager.adata_uuid
        instance_manager_store = self._per_instance_manager_store[self.id]
        instance_manager_store[adata_id] = adata_manager

    def data_registry(self, registry_key: str) -> np.ndarray | pd.DataFrame:
        """Returns the object in AnnData associated with the key in the data registry.

        Parameters
        ----------
        registry_key
            key of object to get from ``self.data_registry``

        Returns
        -------
        The requested data.
        """
        if not self.adata:
            raise ValueError("self.adata is None. Please register AnnData object to access data.")
        else:
            return self._adata_manager.get_from_registry(registry_key)

    def deregister_manager(self, adata: AnnData | None = None):
        """Deregisters the :class:`~scvi.data.AnnDataManager` instance associated with `adata`.

        If `adata` is `None`, deregisters all :class:`~scvi.data.AnnDataManager` instances
        in both the class and instance-specific manager stores, except for the one associated
        with this model instance.
        """
        cls_manager_store = self._setup_adata_manager_store
        instance_manager_store = self._per_instance_manager_store[self.id]

        if adata is None:
            instance_managers_to_clear = list(instance_manager_store.keys())
            cls_managers_to_clear = list(cls_manager_store.keys())
        else:
            adata_manager = self._get_most_recent_anndata_manager(adata, required=True)
            cls_managers_to_clear = [adata_manager.adata_uuid]
            instance_managers_to_clear = [adata_manager.adata_uuid]

        for adata_id in cls_managers_to_clear:
            # don't clear the current manager by default
            is_current_adata = adata is None and adata_id == self.adata_manager.adata_uuid
            if is_current_adata or adata_id not in cls_manager_store:
                continue
            del cls_manager_store[adata_id]

        for adata_id in instance_managers_to_clear:
            # don't clear the current manager by default
            is_current_adata = adata is None and adata_id == self.adata_manager.adata_uuid
            if is_current_adata or adata_id not in instance_manager_store:
                continue
            del instance_manager_store[adata_id]

    @classmethod
    def _get_most_recent_anndata_manager(
        cls, adata: AnnOrMuData, required: bool = False
    ) -> AnnDataManager | None:
        """Retrieves the :class:`~scvi.data.AnnDataManager` for a given AnnData object.

        Checks for the most recent :class:`~scvi.data.AnnDataManager` created for the given AnnData
        object via ``setup_anndata()`` on model initialization. Unlike
        :meth:`scvi.model.base.BaseModelClass.get_anndata_manager`, this method is not model
        instance specific and can be called before a model is fully initialized.

        Parameters
        ----------
        adata
            AnnData object to find manager instance for.
        required
            If True, errors on missing manager. Otherwise, returns None when manager is missing.
        """
        if _SCVI_UUID_KEY not in adata.uns:
            if required:
                raise ValueError(
                    f"Please set up your AnnData with {cls.__name__}.setup_anndata first."
                )
            return None

        adata_id = adata.uns[_SCVI_UUID_KEY]

        if adata_id not in cls._setup_adata_manager_store:
            if required:
                raise ValueError(
                    f"Please set up your AnnData with {cls.__name__}.setup_anndata first. "
                    "It appears the AnnData object has been setup with a different model."
                )
            return None

        adata_manager = cls._setup_adata_manager_store[adata_id]
        if adata_manager.adata is not adata:
            raise ValueError(
                "The provided AnnData object does not match the AnnData object "
                "previously provided for setup. Did you make a copy?"
            )

        return adata_manager

    def get_anndata_manager(
        self, adata: AnnOrMuData, required: bool = False
    ) -> AnnDataManager | None:
        """Retrieves the :class:`~scvi.data.AnnDataManager` for a given AnnData object.

        Requires ``self.id`` has been set. Checks for an :class:`~scvi.data.AnnDataManager`
        specific to this model instance.

        Parameters
        ----------
        adata
            AnnData object to find manager instance for.
        required
            If True, errors on missing manager. Otherwise, returns None when manager is missing.
        """
        cls = self.__class__
        if not adata:
            return None

        if _SCVI_UUID_KEY not in adata.uns:
            if required:
                raise ValueError(
                    f"Please set up your AnnData with {cls.__name__}.setup_anndata'"
                    "or {cls.__name__}.setup_mudata first."
                )
            return None

        adata_id = adata.uns[_SCVI_UUID_KEY]
        if self.id not in cls._per_instance_manager_store:
            if required:
                raise AssertionError(
                    "Unable to find instance specific manager store. "
                    "The model has likely not been initialized with an AnnData object."
                )
            return None
        elif adata_id not in cls._per_instance_manager_store[self.id]:
            if required:
                raise AssertionError(
                    "Please call ``self._validate_anndata`` on this AnnData or MuData object."
                )
            return None

        adata_manager = cls._per_instance_manager_store[self.id][adata_id]
        if adata_manager.adata is not adata:
            logger.info("AnnData object appears to be a copy. Attempting to transfer setup.")
            _assign_adata_uuid(adata, overwrite=True)
            adata_manager = self.adata_manager.transfer_fields(adata)
            self._register_manager_for_instance(adata_manager)

        return adata_manager

    def get_from_registry(
        self,
        adata: AnnOrMuData,
        registry_key: str,
    ) -> np.ndarray:
        """Returns the object in AnnData associated with the key in the data registry.

        AnnData object should be registered with the model prior to calling this function
        via the ``self._validate_anndata`` method.

        Parameters
        ----------
        registry_key
            key of object to get from data registry.
        adata
            AnnData to pull data from.

        Returns
        -------
        The requested data as a NumPy array.
        """
        adata_manager = self.get_anndata_manager(adata)
        if adata_manager is None:
            raise AssertionError(
                "AnnData not registered with model. Call `self._validate_anndata` "
                "prior to calling this function."
            )
        return adata_manager.get_from_registry(registry_key)

    def _make_data_loader(
        self,
        adata: AnnOrMuData,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        shuffle: bool = False,
        data_loader_class=None,
        **data_loader_kwargs,
    ):
        """Create a AnnDataLoader object for data iteration.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        shuffle
            Whether observations are shuffled each iteration though
        data_loader_class
            Class to use for data loader
        data_loader_kwargs
            Kwargs to the class-specific data loader class
        """
        adata_manager = self.get_anndata_manager(adata)
        if adata_manager is None:
            raise AssertionError(
                "AnnDataManager not found. Call `self._validate_anndata` prior to calling this "
                "function."
            )

        adata = adata_manager.adata

        if batch_size is None:
            batch_size = settings.batch_size
        if indices is None:
            indices = np.arange(adata.n_obs)
        if data_loader_class is None:
            data_loader_class = self._data_loader_cls

        if "num_workers" not in data_loader_kwargs:
            data_loader_kwargs.update({"num_workers": settings.dl_num_workers})
        if "persistent_workers" not in data_loader_kwargs:
            data_loader_kwargs.update({"persistent_workers": settings.dl_persistent_workers})

        dl = data_loader_class(
            adata_manager,
            shuffle=shuffle,
            indices=indices,
            batch_size=batch_size,
            **data_loader_kwargs,
        )
        return dl

    def _validate_anndata(
        self, adata: AnnOrMuData | None = None, copy_if_view: bool = True
    ) -> AnnData:
        """Validate anndata has been properly registered, transfer if necessary."""
        if adata is None:
            adata = self.adata

        _check_if_view(adata, copy_if_view=copy_if_view)

        adata_manager = self.get_anndata_manager(adata)
        if adata_manager is None:
            logger.info(
                "Input AnnData not setup with scvi-tools. "
                + "attempting to transfer AnnData setup"
            )
            self._register_manager_for_instance(self.adata_manager.transfer_fields(adata))
        else:
            # Case where correct AnnDataManager is found, replay registration as necessary.
            adata_manager.validate()

        return adata

    def transfer_fields(self, adata: AnnOrMuData, **kwargs) -> AnnData:
        """Transfer fields from a model to an AnnData object."""
        if self.adata:
            return self.adata_manager.transfer_fields(adata, **kwargs)
        else:
            raise ValueError("Model need to be initialized with AnnData to transfer fields.")

    def _check_if_trained(self, warn: bool = True, message: str = _UNTRAINED_WARNING_MESSAGE):
        """Check if the model is trained.

        If not trained and `warn` is True, raise a warning, else raise a RuntimeError.
        """
        if not self.is_trained_:
            if warn:
                warnings.warn(message, UserWarning, stacklevel=settings.warnings_stacklevel)
            else:
                raise RuntimeError(message)

    @property
    def is_trained(self) -> bool:
        """Whether the model has been trained."""
        return self.is_trained_

    @property
    def test_indices(self) -> np.ndarray:
        """Observations that are in test set."""
        return self.test_indices_

    @property
    def train_indices(self) -> np.ndarray:
        """Observations that are in train set."""
        return self.train_indices_

    @property
    def validation_indices(self) -> np.ndarray:
        """Observations that are in validation set."""
        return self.validation_indices_

    @train_indices.setter
    def train_indices(self, value):
        self.train_indices_ = value

    @test_indices.setter
    def test_indices(self, value):
        self.test_indices_ = value

    @validation_indices.setter
    def validation_indices(self, value):
        self.validation_indices_ = value

    @is_trained.setter
    def is_trained(self, value):
        self.is_trained_ = value

    @property
    def history(self):
        """Returns computed metrics during training."""
        return self.history_

    def _get_user_attributes(self):
        """Returns all the self attributes defined in a model class, e.g., `self.is_trained_`."""
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        attributes = [a for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))]
        attributes = [a for a in attributes if not a[0].startswith("_abc_")]
        return attributes

    def _get_init_params(self, locals):
        """Returns the model init signature with associated passed in values.

        Ignores the initial AnnData or Registry.
        """
        init = self.__init__
        sig = inspect.signature(init)
        parameters = sig.parameters.values()

        init_params = [p.name for p in parameters]
        all_params = {p: locals[p] for p in locals if p in init_params}
        all_params = {
            k: v
            for (k, v) in all_params.items()
            if not isinstance(v, AnnData)
            and not isinstance(v, MuData)
            and k not in ("adata", "registry")
        }
        # not very efficient but is explicit
        # separates variable params (**kwargs) from non variable params into two dicts
        non_var_params = [p.name for p in parameters if p.kind != p.VAR_KEYWORD]
        non_var_params = {k: v for (k, v) in all_params.items() if k in non_var_params}
        var_params = [p.name for p in parameters if p.kind == p.VAR_KEYWORD]
        var_params = {k: v for (k, v) in all_params.items() if k in var_params}

        user_params = {"kwargs": var_params, "non_kwargs": non_var_params}

        return user_params

    @abstractmethod
    def train(self):
        """Trains the model."""

    def save(
        self,
        dir_path: str,
        prefix: str | None = None,
        overwrite: bool = False,
        save_anndata: bool = False,
        save_kwargs: dict | None = None,
        legacy_mudata_format: bool = False,
        datamodule: LightningDataModule | None = None,
        **anndata_write_kwargs,
    ):
        """Save the state of the model.

        Neither the trainer optimizer state nor the trainer history are saved.
        Model files are not expected to be reproducibly saved and loaded across versions
        until we reach version 1.0.

        Parameters
        ----------
        dir_path
            Path to a directory.
        prefix
            Prefix to prepend to saved file names.
        overwrite
            Overwrite existing data or not. If `False` and directory
            already exists at `dir_path`, error will be raised.
        save_anndata
            If True, also saves the anndata
        save_kwargs
            Keyword arguments passed into :func:`~torch.save`.
        legacy_mudata_format
            If ``True``, saves the model ``var_names`` in the legacy format if the model was
            trained with a :class:`~mudata.MuData` object. The legacy format is a flat array with
            variable names across all modalities concatenated, while the new format is a dictionary
            with keys corresponding to the modality names and values corresponding to the variable
            names for each modality.
        datamodule
            ``EXPERIMENTAL`` A :class:`~lightning.pytorch.core.LightningDataModule` instance to use
            for training in place of the default :class:`~scvi.dataloaders.DataSplitter`. Can only
            be passed in if the model was not initialized with :class:`~anndata.AnnData`.
        anndata_write_kwargs
            Kwargs for :meth:`~anndata.AnnData.write`
        """
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide another directory for saving."
            )

        file_name_prefix = prefix or ""
        save_kwargs = save_kwargs or {}

        if save_anndata:
            file_suffix = ""
            if isinstance(self.adata, AnnData):
                file_suffix = SAVE_KEYS.ADATA_FNAME
            elif isinstance(self.adata, MuData):
                file_suffix = SAVE_KEYS.MDATA_FNAME
            self.adata.write(
                os.path.join(dir_path, f"{file_name_prefix}{file_suffix}"),
                **anndata_write_kwargs,
            )

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}{SAVE_KEYS.MODEL_FNAME}")

        # save the model state dict and the trainer state dict only
        model_state_dict = self.module.state_dict()
        model_state_dict["pyro_param_store"] = pyro.get_param_store().get_state()

        var_names = self.get_var_names(legacy_mudata_format=legacy_mudata_format)

        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        method_name = self.registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        if method_name == "setup_datamodule":
            user_attributes.update(
                {
                    "n_batch": datamodule.n_batch,
                    "n_extra_categorical_covs": datamodule.registry["field_registries"][
                        "extra_categorical_covs"
                    ]["summary_stats"]["n_extra_categorical_covs"],
                    "n_extra_continuous_covs": datamodule.registry["field_registries"][
                        "extra_continuous_covs"
                    ]["summary_stats"]["n_extra_continuous_covs"],
                    "n_labels": datamodule.n_labels,
                    "n_vars": datamodule.n_vars,
                    "batch_labels": datamodule.batch_labels,
                    "label_keys": datamodule.label_keys,
                }
            )
            if "datamodule" in user_attributes["init_params_"]["non_kwargs"]:
                user_attributes["init_params_"]["non_kwargs"]["datamodule"] = type(
                    user_attributes["init_params_"]["non_kwargs"]["datamodule"]
                ).__name__

        torch.save(
            {
                SAVE_KEYS.MODEL_STATE_DICT_KEY: model_state_dict,
                SAVE_KEYS.VAR_NAMES_KEY: var_names,
                SAVE_KEYS.ATTR_DICT_KEY: user_attributes,
            },
            model_save_path,
            **save_kwargs,
        )

    @classmethod
    @devices_dsp.dedent
    def load(
        cls,
        dir_path: str,
        adata: AnnOrMuData | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        prefix: str | None = None,
        backup_url: str | None = None,
        datamodule: LightningDataModule | None = None,
    ):
        """Instantiate a model from the saved output.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the saved `scvi` setup dictionary.
            If None, will check for and load anndata saved with the model.
            If False, will load the model without AnnData.
        %(param_accelerator)s
        %(param_device)s
        prefix
            Prefix of saved file names.
        backup_url
            URL to retrieve saved outputs from if not present on disk.
        datamodule
            ``EXPERIMENTAL`` A :class:`~lightning.pytorch.core.LightningDataModule` instance to use
            for training in place of the default :class:`~scvi.dataloaders.DataSplitter`. Can only
            be passed in if the model was not initialized with :class:`~anndata.AnnData`.

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> model = ModelClass.load(save_path, adata)
        >>> model.get_....
        """
        load_adata = adata is None
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        (
            attr_dict,
            var_names,
            model_state_dict,
            new_adata,
        ) = _load_saved_files(
            dir_path,
            load_adata,
            map_location=device,
            prefix=prefix,
            backup_url=backup_url,
        )
        adata = new_adata if new_adata is not None else adata

        registry = attr_dict.pop("registry_")
        if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
            raise ValueError("It appears you are loading a model from a different class.")

        # Calling ``setup_anndata`` method with the original arguments passed into
        # the saved model. This enables simple backwards compatibility in the case of
        # newly introduced fields or parameters.
        if adata:
            if _SETUP_ARGS_KEY not in registry:
                raise ValueError(
                    "Saved model does not contain original setup inputs. "
                    "Cannot load the original setup."
                )
            _validate_var_names(adata, var_names)
            method_name = registry.get(_SETUP_METHOD_NAME, "setup_anndata")
            if method_name != "setup_datamodule":
                getattr(cls, method_name)(
                    adata, source_registry=registry, **registry[_SETUP_ARGS_KEY]
                )

        model = _initialize_model(cls, adata, registry, attr_dict, datamodule)
        pyro_param_store = model_state_dict.pop("pyro_param_store", None)

        method_name = registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        if method_name == "setup_datamodule":
            attr_dict["n_input"] = attr_dict["n_vars"]
            module_exp_params = inspect.signature(model._module_cls).parameters.keys()
            common_keys1 = list(attr_dict.keys() & module_exp_params)
            common_keys2 = model.init_params_["non_kwargs"].keys() & module_exp_params
            common_items1 = {key: attr_dict[key] for key in common_keys1}
            common_items2 = {key: model.init_params_["non_kwargs"][key] for key in common_keys2}
            module = model._module_cls(**common_items1, **common_items2)
            model.module = module
        else:
            model.module.on_load(model, pyro_param_store=pyro_param_store)
        model.module.load_state_dict(model_state_dict)

        model.to_device(device)

        model.module.eval()
        if adata:
            model._validate_anndata(adata)
        return model

    @classmethod
    def convert_legacy_save(
        cls,
        dir_path: str,
        output_dir_path: str,
        overwrite: bool = False,
        prefix: str | None = None,
        **save_kwargs,
    ) -> None:
        """Converts a legacy saved model (<v0.15.0) to the updated save format.

        Parameters
        ----------
        dir_path
            Path to directory where legacy model is saved.
        output_dir_path
            Path to save converted save files.
        overwrite
            Overwrite existing data or not. If ``False`` and directory
            already exists at ``output_dir_path``, error will be raised.
        prefix
            Prefix of saved file names.
        **save_kwargs
            Keyword arguments passed into :func:`~torch.save`.
        """
        if not os.path.exists(output_dir_path) or overwrite:
            os.makedirs(output_dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{output_dir_path} already exists. Please provide an unexisting directory for "
                "saving."
            )

        file_name_prefix = prefix or ""
        model_state_dict, var_names, attr_dict, _ = _load_legacy_saved_files(
            dir_path, file_name_prefix, load_adata=False
        )

        if "scvi_setup_dict_" in attr_dict:
            scvi_setup_dict = attr_dict.pop("scvi_setup_dict_")
            unlabeled_category_key = "unlabeled_category_"
            unlabeled_category = attr_dict.get(unlabeled_category_key, None)
            attr_dict["registry_"] = registry_from_setup_dict(
                cls,
                scvi_setup_dict,
                unlabeled_category=unlabeled_category,
            )

        model_save_path = os.path.join(output_dir_path, f"{file_name_prefix}model.pt")
        torch.save(
            {
                SAVE_KEYS.MODEL_STATE_DICT_KEY: model_state_dict,
                SAVE_KEYS.VAR_NAMES_KEY: var_names,
                SAVE_KEYS.ATTR_DICT_KEY: attr_dict,
            },
            model_save_path,
            **save_kwargs,
        )

    @property
    def summary_string(self):
        """Summary string of the model."""
        summary_string = self._model_summary_string
        summary_string += "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        return summary_string

    @property
    def get_normalized_function_name(self) -> str:
        """What the get normalized functions name is"""
        return self.get_normalized_function_name_

    @get_normalized_function_name.setter
    def get_normalized_function_name(self, value):
        self.get_normalized_function_name_ = value

    def __repr__(self):
        rich.print(self.summary_string)
        return ""

    @classmethod
    @abstractmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        *args,
        **kwargs,
    ):
        """%(summary)s.

        Each model class deriving from this class provides parameters to this method
        according to its needs. To operate correctly with the model initialization,
        the implementation must call :meth:`~scvi.model.base.BaseModelClass.register_manager`
        on a model-specific instance of :class:`~scvi.data.AnnDataManager`.
        """

    @staticmethod
    def view_setup_args(dir_path: str, prefix: str | None = None) -> None:
        """Print args used to setup a saved model.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        prefix
            Prefix of saved file names.
        """
        registry = BaseModelClass.load_registry(dir_path, prefix)
        AnnDataManager.view_setup_method_args(registry)

    @staticmethod
    def load_registry(dir_path: str, prefix: str | None = None) -> dict:
        """Return the full registry saved with the model.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        prefix
            Prefix of saved file names.

        Returns
        -------
        The full registry saved with the model
        """
        attr_dict = _load_saved_files(dir_path, False, prefix=prefix)[0]

        # Legacy support for old setup dict format.
        if "scvi_setup_dict_" in attr_dict:
            raise NotImplementedError(
                "Viewing setup args for pre v0.15.0 models is unsupported. "
                "Update your save files with ``convert_legacy_save`` first."
            )

        return attr_dict.pop("registry_")

    def view_anndata_setup(
        self, adata: AnnOrMuData | None = None, hide_state_registries: bool = False
    ) -> None:
        """Print summary of the setup for the initial AnnData or a given AnnData object.

        Parameters
        ----------
        adata
            AnnData object setup with ``setup_anndata`` or
            :meth:`~scvi.data.AnnDataManager.transfer_fields`.
        hide_state_registries
            If True, prints a shortened summary without details of each state registry.
        """
        if adata is None:
            adata = self.adata
        try:
            adata_manager = self.get_anndata_manager(adata, required=True)
        except ValueError as err:
            raise ValueError(
                f"Given AnnData not setup with {self.__class__.__name__}. "
                "Cannot view setup summary."
            ) from err
        adata_manager.view_registry(hide_state_registries=hide_state_registries)

    def view_setup_method_args(self) -> None:
        """Prints setup kwargs used to produce a given registry.

        Parameters
        ----------
        registry
            Registry produced by an AnnDataManager.
        """
        model_name = self.registry_[_MODEL_NAME_KEY]
        setup_args = self.registry_[_SETUP_ARGS_KEY]
        if model_name is not None and setup_args is not None:
            rich.print(f"Setup via `{model_name}.setup_anndata` with arguments:")
            rich.pretty.pprint(setup_args)
            rich.print()

    def view_registry(self, hide_state_registries: bool = False) -> None:
        """Prints summary of the registry.

        Parameters
        ----------
        hide_state_registries
            If True, prints a shortened summary without details of each state registry.
        """
        version = self.registry_[_SCVI_VERSION_KEY]
        rich.print(f"Anndata setup with scvi-tools version {version}.")
        rich.print()
        self.view_setup_method_args(self._registry)

        in_colab = "google.colab" in sys.modules
        force_jupyter = None if not in_colab else True
        console = rich.console.Console(force_jupyter=force_jupyter)

        ss = AnnDataManager._get_summary_stats_from_registry(self._registry)
        dr = self._get_data_registry_from_registry(self._registry)
        console.print(self._view_summary_stats(ss))
        console.print(self._view_data_registry(dr))

        if not hide_state_registries:
            for field in self.fields:
                state_registry = self.get_state_registry(field.registry_key)
                t = field.view_state_registry(state_registry)
                if t is not None:
                    console.print(t)

    def get_state_registry(self, registry_key: str) -> attrdict:
        """Returns the state registry for the AnnDataField registered with this instance."""
        return attrdict(self.registry_[_FIELD_REGISTRIES_KEY][registry_key][_STATE_REGISTRY_KEY])

    def get_setup_arg(self, setup_arg: str) -> attrdict:
        """Returns the string provided to setup of a specific setup_arg."""
        return self.registry_[_SETUP_ARGS_KEY][setup_arg]

    @staticmethod
    def _view_summary_stats(
        summary_stats: attrdict, as_markdown: bool = False
    ) -> rich.table.Table | str:
        """Prints summary stats."""
        if not as_markdown:
            t = rich.table.Table(title="Summary Statistics")
        else:
            t = rich.table.Table(box=box.MARKDOWN)

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
        for stat_key, count in summary_stats.items():
            t.add_row(stat_key, str(count))

        if as_markdown:
            console = Console(file=StringIO(), force_jupyter=False)
            console.print(t)
            return console.file.getvalue().strip()

        return t

    @staticmethod
    def _view_data_registry(
        data_registry: attrdict, as_markdown: bool = False
    ) -> rich.table.Table | str:
        """Prints data registry."""
        if not as_markdown:
            t = rich.table.Table(title="Data Registry")
        else:
            t = rich.table.Table(box=box.MARKDOWN)

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

        for registry_key, data_loc in data_registry.items():
            mod_key = getattr(data_loc, _constants._DR_MOD_KEY, None)
            attr_name = data_loc.attr_name
            attr_key = data_loc.attr_key
            scvi_data_str = "adata"
            if mod_key is not None:
                scvi_data_str += f".mod['{mod_key}']"
            if attr_key is None:
                scvi_data_str += f".{attr_name}"
            else:
                scvi_data_str += f".{attr_name}['{attr_key}']"
            t.add_row(registry_key, scvi_data_str)

        if as_markdown:
            console = Console(file=StringIO(), force_jupyter=False)
            console.print(t)
            return console.file.getvalue().strip()

        return t

    def update_setup_method_args(self, setup_method_args: dict):
        """Update setup method args.

        Parameters
        ----------
        setup_method_args
            This is a bit of a misnomer, this is a dict representing kwargs
            of the setup method that will be used to update the existing values
            in the registry of this instance.
        """
        self._registry[_SETUP_ARGS_KEY].update(setup_method_args)

    def get_normalized_expression(self, *args, **kwargs):
        msg = f"get_normalized_expression is not implemented for {self.__class__.__name__}."
        raise NotImplementedError(msg)


class BaseMinifiedModeModelClass(BaseModelClass):
    """Abstract base class for scvi-tools models that can handle minified data."""

    @property
    def minified_data_type(self) -> MinifiedDataType | None:
        """The type of minified data associated with this model, if applicable."""
        if self.adata_manager:
            return (
                self.adata_manager.get_from_registry(REGISTRY_KEYS.MINIFY_TYPE_KEY)
                if REGISTRY_KEYS.MINIFY_TYPE_KEY in self.adata_manager.data_registry
                else None
            )
        else:
            return None

    def minify_adata(
        self,
        minified_data_type: MinifiedDataType = ADATA_MINIFY_TYPE.LATENT_POSTERIOR,
        use_latent_qzm_key: str = "X_latent_qzm",
        use_latent_qzv_key: str = "X_latent_qzv",
    ) -> None:
        """Minify the model's :attr:`~scvi.model.base.BaseModelClass.adata`.

        Minifies the :class:`~anndata.AnnData` object associated with the model according to the
        method specified by ``minified_data_type`` and registers the new fields with the model's
        :class:`~scvi.data.AnnDataManager`. This also sets the ``minified_data_type`` attribute
        of the underlying :class:`~scvi.module.base.BaseModuleClass` instance.

        Parameters
        ----------
        minified_data_type
            Method for minifying the data. One of the following:

            - ``"latent_posterior_parameters"``: Store the latent posterior mean and variance in
                :attr:`~anndata.AnnData.obsm` using the keys ``use_latent_qzm_key`` and
                ``use_latent_qzv_key``.
            - ``"latent_posterior_parameters_with_counts"``: Store the latent posterior mean and
                variance in :attr:`~anndata.AnnData.obsm` using the keys ``use_latent_qzm_key`` and
                ``use_latent_qzv_key``, and the raw count data in :attr:`~anndata.AnnData.X`.
        use_latent_qzm_key
            Key to use for storing the latent posterior mean in :attr:`~anndata.AnnData.obsm` when
            ``minified_data_type`` is ``"latent_posterior"``.
        use_latent_qzv_key
            Key to use for storing the latent posterior variance in :attr:`~anndata.AnnData.obsm`
            when ``minified_data_type`` is ``"latent_posterior"``.

        Notes
        -----
        The modification is not done inplace -- instead the model is assigned a new (minified)
        version of the :class:`~anndata.AnnData`.
        """
        if minified_data_type not in ADATA_MINIFY_TYPE:
            raise NotImplementedError(
                f"Minification method {minified_data_type} is not supported."
            )
        elif not getattr(self.module, "use_observed_lib_size", True):
            raise ValueError(
                "Minification is not supported for models that do not use observed library size."
            )

        keep_count_data = minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR_WITH_COUNTS
        mini_adata = get_minified_adata_scrna(
            adata_manager=self.adata_manager,
            keep_count_data=keep_count_data,
        )
        del mini_adata.uns[_SCVI_UUID_KEY]
        mini_adata.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = minified_data_type
        mini_adata.obsm[self._LATENT_QZM_KEY] = self.adata.obsm[use_latent_qzm_key]
        mini_adata.obsm[self._LATENT_QZV_KEY] = self.adata.obsm[use_latent_qzv_key]
        mini_adata.obs[self._OBSERVED_LIB_SIZE_KEY] = np.squeeze(
            np.asarray(self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY).sum(axis=-1))
        )
        self._update_adata_and_manager_post_minification(
            mini_adata,
            minified_data_type,
        )
        self.module.minified_data_type = minified_data_type

    @classmethod
    def _get_fields_for_adata_minification(
        cls,
        minified_data_type: MinifiedDataType,
    ):
        """Return the fields required for minification of the given type."""
        if minified_data_type not in ADATA_MINIFY_TYPE:
            raise NotImplementedError(
                f"Minification method {minified_data_type} is not supported."
            )

        mini_fields = [
            fields.ObsmField(REGISTRY_KEYS.LATENT_QZM_KEY, cls._LATENT_QZM_KEY),
            fields.ObsmField(REGISTRY_KEYS.LATENT_QZV_KEY, cls._LATENT_QZV_KEY),
            fields.NumericalObsField(REGISTRY_KEYS.OBSERVED_LIB_SIZE, cls._OBSERVED_LIB_SIZE_KEY),
            fields.StringUnsField(REGISTRY_KEYS.MINIFY_TYPE_KEY, _ADATA_MINIFY_TYPE_UNS_KEY),
        ]

        return mini_fields

    def _update_adata_and_manager_post_minification(
        self,
        minified_adata: AnnOrMuData,
        minified_data_type: MinifiedDataType,
    ):
        """Update the :class:`~anndata.AnnData` and :class:`~scvi.data.AnnDataManager` in-place.

        Parameters
        ----------
        minified_adata
            Minified version of :attr:`~scvi.model.base.BaseModelClass.adata`.
        minified_data_type
            Method used for minifying the data.
        keep_count_data
            If ``True``, the full count matrix is kept in the minified
            :attr:`~scvi.model.base.BaseModelClass.adata`.
        """
        self._validate_anndata(minified_adata)
        new_adata_manager = self.get_anndata_manager(minified_adata, required=True)
        new_adata_manager.register_new_fields(
            self._get_fields_for_adata_minification(minified_data_type)
        )
        self.adata = minified_adata

    @property
    def summary_string(self):
        """Summary string of the model."""
        summary_string = super().summary_string
        summary_string += "\nModel's adata is minified?: {}".format(
            hasattr(self, "minified_data_type") and self.minified_data_type is not None
        )
        return summary_string


class BaseMudataMinifiedModeModelClass(BaseModelClass):
    """Abstract base class for scvi-tools models that can handle minified data."""

    @property
    def minified_data_type(self) -> MinifiedDataType | None:
        """The type of minified data associated with this model, if applicable."""
        return (
            self.adata_manager.get_from_registry(REGISTRY_KEYS.MINIFY_TYPE_KEY)
            if REGISTRY_KEYS.MINIFY_TYPE_KEY in self.adata_manager.data_registry
            else None
        )

    def minify_mudata(
        self,
        minified_data_type: MinifiedDataType = ADATA_MINIFY_TYPE.LATENT_POSTERIOR,
        use_latent_qzm_key: str = "X_latent_qzm",
        use_latent_qzv_key: str = "X_latent_qzv",
    ) -> None:
        """Minify the model's :attr:`~scvi.model.base.BaseModelClass.adata`.

        Minifies the :class:`~mudata.MuData` object associated with the model according to the
        method specified by ``minified_data_type`` and registers the new fields with the model's
        :class:`~scvi.data.AnnDataManager`. This also sets the ``minified_data_type`` attribute
        of the underlying :class:`~scvi.module.base.BaseModuleClass` instance.

        Parameters
        ----------
        minified_data_type
            Method for minifying the data. One of the following:

            - ``"latent_posterior_parameters"``: Store the latent posterior mean and variance in
                :attr:`~mudata.MuData.obsm` using the keys ``use_latent_qzm_key`` and
                ``use_latent_qzv_key``.
            - ``"latent_posterior_parameters_with_counts"``: Store the latent posterior mean and
                variance in :attr:`~mudata.MuData.obsm` using the keys ``use_latent_qzm_key`` and
                ``use_latent_qzv_key``, and the raw count data in :attr:`~mudata.MuData.X`.
        use_latent_qzm_key
            Key to use for storing the latent posterior mean in :attr:`~mudata.MuData.obsm` when
            ``minified_data_type`` is ``"latent_posterior"``.
        use_latent_qzv_key
            Key to use for storing the latent posterior variance in :attr:`~mudata.MuData.obsm`
            when ``minified_data_type`` is ``"latent_posterior"``.

        Notes
        -----
        The modification is not done inplace -- instead the model is assigned a new (minified)
        version of the :class:`~mudata.MuData`.
        """
        if self.adata_manager._registry["setup_method_name"] != "setup_mudata":
            raise ValueError(
                f"MuData must be registered with {self.__name__}.setup_mudata to use this method."
            )
        if minified_data_type not in ADATA_MINIFY_TYPE:
            raise NotImplementedError(
                f"Minification method {minified_data_type} is not supported."
            )
        elif not getattr(self.module, "use_observed_lib_size", True):
            raise ValueError(
                "Minification is not supported for models that do not use observed library size."
            )

        keep_count_data = minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR_WITH_COUNTS
        mini_adata = get_minified_mudata(
            adata_manager=self.adata_manager,
            keep_count_data=keep_count_data,
        )
        del mini_adata.uns[_SCVI_UUID_KEY]
        mini_adata.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = minified_data_type
        mini_adata.obsm[self._LATENT_QZM_KEY] = self.adata.obsm[use_latent_qzm_key]
        mini_adata.obsm[self._LATENT_QZV_KEY] = self.adata.obsm[use_latent_qzv_key]
        mini_adata.obs[self._OBSERVED_LIB_SIZE_KEY] = np.squeeze(
            np.asarray(self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY).sum(axis=-1))
        )
        self._update_mudata_and_manager_post_minification(
            mini_adata,
            minified_data_type,
        )
        self.module.minified_data_type = minified_data_type

    @classmethod
    def _get_fields_for_mudata_minification(
        cls,
        minified_data_type: MinifiedDataType,
    ):
        """Return the fields required for minification of the given type."""
        if minified_data_type not in ADATA_MINIFY_TYPE:
            raise NotImplementedError(
                f"Minification method {minified_data_type} is not supported."
            )

        mini_fields = [
            fields.ObsmField(REGISTRY_KEYS.LATENT_QZM_KEY, cls._LATENT_QZM_KEY),
            fields.ObsmField(REGISTRY_KEYS.LATENT_QZV_KEY, cls._LATENT_QZV_KEY),
            fields.NumericalObsField(REGISTRY_KEYS.OBSERVED_LIB_SIZE, cls._OBSERVED_LIB_SIZE_KEY),
            fields.StringUnsField(REGISTRY_KEYS.MINIFY_TYPE_KEY, _ADATA_MINIFY_TYPE_UNS_KEY),
        ]

        return mini_fields

    def _update_mudata_and_manager_post_minification(
        self, minified_adata: AnnOrMuData, minified_data_type: MinifiedDataType
    ):
        """Update the mudata and manager inplace after creating a minified adata."""
        # Register this new adata with the model, creating a new manager in the cache
        self._validate_anndata(minified_adata)
        new_adata_manager = self.get_anndata_manager(minified_adata, required=True)
        # This inplace edits the manager
        new_adata_manager.register_new_fields(
            self._get_fields_for_mudata_minification(minified_data_type),
        )
        new_adata_manager.registry["setup_method_name"] = "setup_mudata"
        # We set the adata attribute of the model as this will update self.registry_
        # and self.adata_manager with the new adata manager
        self.adata = minified_adata

    @property
    def summary_string(self):
        """Summary string of the model."""
        summary_string = super().summary_string
        summary_string += "\nModel's adata is minified?: {}".format(
            hasattr(self, "minified_data_type") and self.minified_data_type is not None
        )
        return summary_string
