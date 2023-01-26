import inspect
import logging
import os
import warnings
from abc import ABCMeta, abstractmethod
from typing import Dict, Optional, Sequence, Type, Union
from uuid import uuid4

import numpy as np
import rich
import torch
from anndata import AnnData
from mudata import MuData

from scvi import REGISTRY_KEYS, settings
from scvi._types import AnnOrMuData, LatentDataType
from scvi.autotune._types import TunableMixin
from scvi.data import AnnDataManager
from scvi.data._compat import registry_from_setup_dict
from scvi.data._constants import (
    _MODEL_NAME_KEY,
    _SCVI_UUID_KEY,
    _SETUP_ARGS_KEY,
    _SETUP_METHOD_NAME,
)
from scvi.data._utils import _assign_adata_uuid, _check_if_view, _get_latent_adata_type
from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_use_gpu_arg
from scvi.model.base._utils import _load_legacy_saved_files
from scvi.utils import attrdict, setup_anndata_dsp

from ._utils import _initialize_model, _load_saved_files, _validate_var_names

logger = logging.getLogger(__name__)


_UNTRAINED_WARNING_MESSAGE = "Trying to query inferred values from an untrained model. Please train the model first."

_SETUP_INPUTS_EXCLUDED_PARAMS = {"adata", "mdata", "kwargs"}


class BaseModelMetaClass(ABCMeta):
    """
    Metaclass for :class:`~scvi.model.base.BaseModelClass`.

    Constructs model class-specific mappings for :class:`~scvi.data.AnnDataManager` instances.
    ``cls._setup_adata_manager_store`` maps from AnnData object UUIDs to :class:`~scvi.data.AnnDataManager` instances.
    This mapping is populated everytime ``cls.setup_anndata()`` is called.
    ``cls._per_isntance_manager_store`` maps from model instance UUIDs to AnnData UUID::class:`~scvi.data.AnnDataManager` mappings.
    These :class:`~scvi.data.AnnDataManager` instances are tied to a single model instance and populated either
    during model initialization or after running ``self._validate_anndata()``.
    """

    @abstractmethod
    def __init__(cls, name, bases, dct):
        cls._setup_adata_manager_store: Dict[
            str, Type[AnnDataManager]
        ] = dict()  # Maps adata id to AnnDataManager instances.
        cls._per_instance_manager_store: Dict[
            str, Dict[str, Type[AnnDataManager]]
        ] = dict()  # Maps model instance id to AnnDataManager mappings.
        super().__init__(name, bases, dct)


class BaseModelClass(TunableMixin, metaclass=BaseModelMetaClass):
    """Abstract class for scvi-tools models."""

    def __init__(self, adata: Optional[AnnOrMuData] = None):
        # check if the given adata is in latent mode and check if the model being created
        # supports latent mode (i.e. inherits from the abstract BaseLatentModeModelClass).
        # If not, raise an error to inform the user of the lack of latent mode functionality
        # for this model
        latent_adata = adata is not None and _get_latent_adata_type(adata) is not None
        if latent_adata and not issubclass(type(self), BaseLatentModeModelClass):
            raise NotImplementedError(
                f"Latent mode currently not supported for the {type(self).__name__} model."
            )
        self.id = str(uuid4())  # Used for cls._manager_store keys.
        if adata is not None:
            self._adata = adata
            self._adata_manager = self._get_most_recent_anndata_manager(
                adata, required=True
            )
            self._register_manager_for_instance(self.adata_manager)
            # Suffix registry instance variable with _ to include it when saving the model.
            self.registry_ = self._adata_manager.registry
            self.summary_stats = self._adata_manager.summary_stats

        self.is_trained_ = False
        self._model_summary_string = ""
        self.train_indices_ = None
        self.test_indices_ = None
        self.validation_indices_ = None
        self.history_ = None
        self._data_loader_cls = AnnDataLoader

    @property
    def adata(self) -> AnnOrMuData:
        """Data attached to model instance."""
        return self._adata

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

    def to_device(self, device: Union[str, int]):
        """
        Move model to device.

        Parameters
        ----------
        device
            Device to move model to. Options: 'cpu' for CPU, integer GPU index (eg. 0),
            or 'cuda:X' where X is the GPU index (eg. 'cuda:0'). See torch.device for more info.

        Examples
        --------
        >>> adata = scvi.data.synthetic_iid()
        >>> model = scvi.model.SCVI(adata)
        >>> model.to_device('cpu')      # moves model to CPU
        >>> model.to_device('cuda:0')   # moves model to GPU 0
        >>> model.to_device(0)          # also moves model to GPU 0
        """
        my_device = torch.device(device)
        self.module.to(my_device)

    @property
    def device(self) -> str:
        """The current device that the module's params are on."""
        return self.module.device

    @staticmethod
    def _get_setup_method_args(**setup_locals) -> dict:
        """
        Returns a dictionary organizing the arguments used to call ``setup_anndata``.

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
        setup_args = dict()
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
        modalities: Dict[str, str], setup_method_args: dict
    ) -> attrdict:
        """
        Preprocesses a ``modalities`` dictionary used in ``setup_mudata()`` to map modality names.

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
            raise ValueError(
                f"Extraneous modality mapping(s) detected: {extra_modalities}"
            )
        return attrdict(filtered_modalities)

    @classmethod
    def register_manager(cls, adata_manager: AnnDataManager):
        """
        Registers an :class:`~scvi.data.AnnDataManager` instance with this model class.

        Stores the :class:`~scvi.data.AnnDataManager` reference in a class-specific manager store.
        Intended for use in the ``setup_anndata()`` class method followed up by retrieval of the
        :class:`~scvi.data.AnnDataManager` via the ``_get_most_recent_anndata_manager()`` method in
        the model init method.

        Notes
        -----
        Subsequent calls to this method with an :class:`~scvi.data.AnnDataManager` instance referring to the same
        underlying AnnData object will overwrite the reference to previous :class:`~scvi.data.AnnDataManager`.
        """
        adata_id = adata_manager.adata_uuid
        cls._setup_adata_manager_store[adata_id] = adata_manager

    def _register_manager_for_instance(self, adata_manager: AnnDataManager):
        """
        Registers an :class:`~scvi.data.AnnDataManager` instance with this model instance.

        Creates a model-instance specific mapping in ``cls._per_instance_manager_store`` for this
        :class:`~scvi.data.AnnDataManager` instance.
        """
        if self.id not in self._per_instance_manager_store:
            self._per_instance_manager_store[self.id] = dict()

        adata_id = adata_manager.adata_uuid
        instance_manager_store = self._per_instance_manager_store[self.id]
        instance_manager_store[adata_id] = adata_manager

    @classmethod
    def _get_most_recent_anndata_manager(
        cls, adata: AnnOrMuData, required: bool = False
    ) -> Optional[AnnDataManager]:
        """
        Retrieves the :class:`~scvi.data.AnnDataManager` for a given AnnData object specific to this model class.

        Checks for the most recent :class:`~scvi.data.AnnDataManager` created for the given AnnData object via
        ``setup_anndata()`` on model initialization. Unlike :meth:`scvi.model.base.BaseModelClass.get_anndata_manager`,
        this method is not model instance specific and can be called before a model is fully initialized.

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
    ) -> Optional[AnnDataManager]:
        """
        Retrieves the :class:`~scvi.data.AnnDataManager` for a given AnnData object specific to this model instance.

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
        if _SCVI_UUID_KEY not in adata.uns:
            if required:
                raise ValueError(
                    f"Please set up your AnnData with {cls.__name__}.setup_anndata first."
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
                    "Please call ``self._validate_anndata`` on this AnnData object."
                )
            return None

        adata_manager = cls._per_instance_manager_store[self.id][adata_id]
        if adata_manager.adata is not adata:
            logger.info(
                "AnnData object appears to be a copy. Attempting to transfer setup."
            )
            _assign_adata_uuid(adata, overwrite=True)
            adata_manager = self.adata_manager.transfer_fields(adata)
            self._register_manager_for_instance(adata_manager)

        return adata_manager

    def get_from_registry(
        self,
        adata: AnnOrMuData,
        registry_key: str,
    ) -> np.ndarray:
        """
        Returns the object in AnnData associated with the key in the data registry.

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
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        shuffle: bool = False,
        data_loader_class=None,
        **data_loader_kwargs,
    ):
        """
        Create a AnnDataLoader object for data iteration.

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
                "AnnDataManager not found. Call `self._validate_anndata` prior to calling this function."
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

        dl = data_loader_class(
            adata_manager,
            shuffle=shuffle,
            indices=indices,
            batch_size=batch_size,
            **data_loader_kwargs,
        )
        return dl

    def _validate_anndata(
        self, adata: Optional[AnnOrMuData] = None, copy_if_view: bool = True
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
            self._register_manager_for_instance(
                self.adata_manager.transfer_fields(adata)
            )
        else:
            # Case where correct AnnDataManager is found, replay registration as necessary.
            adata_manager.validate()

        return adata

    def _check_if_trained(
        self, warn: bool = True, message: str = _UNTRAINED_WARNING_MESSAGE
    ):
        """
        Check if the model is trained.

        If not trained and `warn` is True, raise a warning, else raise a RuntimeError.
        """
        if not self.is_trained_:
            if warn:
                warnings.warn(message)
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
        attributes = [
            a for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))
        ]
        attributes = [a for a in attributes if not a[0].startswith("_abc_")]
        return attributes

    def _get_init_params(self, locals):
        """
        Returns the model init signature with associated passed in values.

        Ignores the initial AnnData.
        """
        init = self.__init__
        sig = inspect.signature(init)
        parameters = sig.parameters.values()

        init_params = [p.name for p in parameters]
        all_params = {p: locals[p] for p in locals if p in init_params}
        all_params = {
            k: v
            for (k, v) in all_params.items()
            if not isinstance(v, AnnData) and not isinstance(v, MuData)
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
        prefix: Optional[str] = None,
        overwrite: bool = False,
        save_anndata: bool = False,
        **anndata_write_kwargs,
    ):
        """
        Save the state of the model.

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
        anndata_write_kwargs
            Kwargs for :meth:`~anndata.AnnData.write`
        """
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide another directory for saving.".format(
                    dir_path
                )
            )

        file_name_prefix = prefix or ""
        if save_anndata:
            file_suffix = ""
            if isinstance(self.adata, AnnData):
                file_suffix = "adata.h5ad"
            elif isinstance(self.adata, MuData):
                file_suffix = "mdata.h5mu"
            self.adata.write(
                os.path.join(dir_path, f"{file_name_prefix}{file_suffix}"),
                **anndata_write_kwargs,
            )

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")

        # save the model state dict and the trainer state dict only
        model_state_dict = self.module.state_dict()

        var_names = self.adata.var_names.astype(str)
        var_names = var_names.to_numpy()

        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        torch.save(
            dict(
                model_state_dict=model_state_dict,
                var_names=var_names,
                attr_dict=user_attributes,
            ),
            model_save_path,
        )

    @classmethod
    def load(
        cls,
        dir_path: str,
        adata: Optional[AnnOrMuData] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        prefix: Optional[str] = None,
        backup_url: Optional[str] = None,
    ):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the saved `scvi` setup dictionary.
            If None, will check for and load anndata saved with the model.
        use_gpu
            Load model on default GPU if available (if None or True),
            or index of GPU to use (if int), or name of GPU (if str), or use CPU (if False).
        prefix
            Prefix of saved file names.
        backup_url
            URL to retrieve saved outputs from if not present on disk.

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> model = ModelClass.load(save_path, adata) # use the name of the model class used to save
        >>> model.get_....
        """
        load_adata = adata is None
        _, _, device = parse_use_gpu_arg(use_gpu)

        (attr_dict, var_names, model_state_dict, new_adata,) = _load_saved_files(
            dir_path,
            load_adata,
            map_location=device,
            prefix=prefix,
            backup_url=backup_url,
        )
        adata = new_adata if new_adata is not None else adata

        _validate_var_names(adata, var_names)

        registry = attr_dict.pop("registry_")
        if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
            raise ValueError(
                "It appears you are loading a model from a different class."
            )

        if _SETUP_ARGS_KEY not in registry:
            raise ValueError(
                "Saved model does not contain original setup inputs. "
                "Cannot load the original setup."
            )

        # Calling ``setup_anndata`` method with the original arguments passed into
        # the saved model. This enables simple backwards compatibility in the case of
        # newly introduced fields or parameters.
        method_name = registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        getattr(cls, method_name)(
            adata, source_registry=registry, **registry[_SETUP_ARGS_KEY]
        )

        model = _initialize_model(cls, adata, attr_dict)
        model.module.on_load(model)
        model.module.load_state_dict(model_state_dict)

        model.to_device(device)
        model.module.eval()
        model._validate_anndata(adata)
        return model

    @classmethod
    def convert_legacy_save(
        cls,
        dir_path: str,
        output_dir_path: str,
        overwrite: bool = False,
        prefix: Optional[str] = None,
    ) -> None:
        """
        Converts a legacy saved model (<v0.15.0) to the updated save format.

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
        """
        if not os.path.exists(output_dir_path) or overwrite:
            os.makedirs(output_dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{output_dir_path} already exists. Please provide an unexisting directory for saving."
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
            dict(
                model_state_dict=model_state_dict,
                var_names=var_names,
                attr_dict=attr_dict,
            ),
            model_save_path,
        )

    def __repr__(
        self,
    ):
        summary_string = self._model_summary_string
        summary_string += "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        summary_string += "\nLatent model? {}".format(
            hasattr(self, "latent_data_type") and self.latent_data_type is not None
        )
        rich.print(summary_string)

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
        """
        %(summary)s.

        Each model class deriving from this class provides parameters to this method
        according to its needs. To operate correctly with the model initialization,
        the implementation must call :meth:`~scvi.model.base.BaseModelClass.register_manager`
        on a model-specific instance of :class:`~scvi.data.AnnDataManager`.
        """

    @staticmethod
    def view_setup_args(dir_path: str, prefix: Optional[str] = None) -> None:
        """
        Print args used to setup a saved model.

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
    def load_registry(dir_path: str, prefix: Optional[str] = None) -> dict:
        """
        Return the full registry saved with the model.

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
        self, adata: Optional[AnnOrMuData] = None, hide_state_registries: bool = False
    ) -> None:
        """
        Print summary of the setup for the initial AnnData or a given AnnData object.

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
        except ValueError:
            raise ValueError(
                f"Given AnnData not setup with {self.__class__.__name__}. "
                "Cannot view setup summary."
            )
        adata_manager.view_registry(hide_state_registries=hide_state_registries)


class BaseLatentModeModelClass(BaseModelClass):
    """Abstract base class for scvi-tools models that support latent mode."""

    @property
    def latent_data_type(self) -> Union[LatentDataType, None]:
        """The latent data type associated with this model."""
        return (
            self.adata_manager.get_from_registry(REGISTRY_KEYS.LATENT_MODE_KEY)
            if REGISTRY_KEYS.LATENT_MODE_KEY in self.adata_manager.data_registry
            else None
        )

    @abstractmethod
    def to_latent_mode(
        self,
        *args,
        **kwargs,
    ):
        """
        Put the model into latent mode.

        The model is put into latent mode by registering new anndata fields
        required for latent mode support (can be model-specific) and marking
        the module as latent. Note that this modifies the anndata (and subsequently
        the model and module properties) in place. Please make a copy of those objects
        (before calling this function) if needed.
        """

    @staticmethod
    @abstractmethod
    def _get_latent_fields(mode: LatentDataType):
        """Return the anndata fields required for latent mode support."""
