import torch
import pickle
import os
import logging
import inspect
import numpy as np

from anndata import AnnData
from scvi.data._utils import (
    _check_anndata_setup_equivalence,
    _check_nonnegative_integers,
)
from scvi import _CONSTANTS, settings
from typing import Optional, Sequence
from abc import ABC, abstractmethod
from scvi.data import get_from_registry, transfer_anndata_setup

logger = logging.getLogger(__name__)


class BaseModelClass(ABC):
    def __init__(self, adata: Optional[AnnData] = None, use_cuda=False):
        if adata is not None:
            if "_scvi" not in adata.uns.keys():
                raise ValueError(
                    "Please setup your AnnData with scvi.data.setup_anndata(adata) first"
                )
            self.adata = adata
            self.scvi_setup_dict_ = adata.uns["_scvi"]
            self.summary_stats = self.scvi_setup_dict_["summary_stats"]
            self._validate_anndata(adata, copy_if_view=False)

        self.is_trained_ = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self._model_summary_string = ""
        self.train_indices_ = None
        self.test_indices_ = None
        self.validation_indices_ = None

    def _make_scvi_dl(
        self,
        adata: AnnData,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        **data_loader_kwargs,
    ):
        """Create a ScviDataLoader object for data iteration."""
        if batch_size is None:
            batch_size = settings.batch_size
        if indices is None:
            indices = np.arange(adata.n_obs)
        post = self._scvi_dl_class(
            self.model,
            adata,
            shuffle=False,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=batch_size,
            **data_loader_kwargs,
        ).sequential()
        return post

    def _validate_anndata(
        self, adata: Optional[AnnData] = None, copy_if_view: bool = True
    ):
        """Validate anndata has been properly registered, transfer if necessary."""
        if adata is None:
            adata = self.adata
        if adata.is_view:
            if copy_if_view:
                logger.info("Received view of anndata, making copy.")
                adata = adata.copy()
            else:
                raise ValueError("Please run `adata = adata.copy()`")

        if "_scvi" not in adata.uns_keys():
            logger.info(
                "Input adata not setup with scvi. "
                + "attempting to transfer anndata setup"
            )
            transfer_anndata_setup(self.scvi_setup_dict_, adata)
        is_nonneg_int = _check_nonnegative_integers(
            get_from_registry(adata, _CONSTANTS.X_KEY)
        )
        if not is_nonneg_int:
            logger.warning(
                "Make sure the registered X field in anndata contains unnormalized count data."
            )

        _check_anndata_setup_equivalence(self.scvi_setup_dict_, adata)

        return adata

    @property
    @abstractmethod
    def _scvi_dl_class(self):
        pass

    @property
    @abstractmethod
    def _trainer_class(self):
        pass

    @abstractmethod
    def train(self):
        pass

    @property
    def is_trained(self):
        return self.is_trained_

    @property
    def test_indices(self):
        return self.test_indices_

    @property
    def train_indices(self):
        return self.train_indices_

    @property
    def validation_indices(self):
        return self.validation_indices_

    def _get_user_attributes(self):
        # returns all the self attributes defined in a model class, eg, self.is_trained_
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        attributes = [
            a for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))
        ]
        attributes = [a for a in attributes if not a[0].startswith("_abc_")]
        return attributes

    def _get_init_params(self, locals):
        # returns the model init signiture with associated passed in values
        # except the anndata objects passed in
        init = self.__init__
        sig = inspect.signature(init)
        init_params = [p for p in sig.parameters]
        user_params = {p: locals[p] for p in locals if p in init_params}
        user_params = {
            k: v for (k, v) in user_params.items() if not isinstance(v, AnnData)
        }
        return user_params

    def save(self, dir_path: str, overwrite: bool = False):
        """
        Save the state of the model.

        Neither the trainer optimizer state nor the trainer history are saved.
        Model files are not expected to be reproducibly saved and loaded across versions
        until we reach version 1.0.

        Parameters
        ----------
        dir_path
            Path to a directory.
        overwrite
            Overwrite existing data or not. If `False` and directory
            already exists at `dir_path`, error will be raised.
        """
        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}
        # save the model state dict and the trainer state dict only
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )
        torch.save(self.model.state_dict(), os.path.join(dir_path, "model_params.pt"))
        with open(os.path.join(dir_path, "attr.pkl"), "wb") as f:
            pickle.dump(user_attributes, f)

    @classmethod
    def load(cls, adata: AnnData, dir_path: str, use_cuda: bool = False):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :func:`~scvi.data.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
        dir_path
            Path to saved outputs.
        use_cuda
            Whether to load model on GPU.

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> vae = SCVI.load(adata, save_path)
        >>> vae.get_latent_representation()
        """
        model_path = os.path.join(dir_path, "model_params.pt")
        setup_dict_path = os.path.join(dir_path, "attr.pkl")
        with open(setup_dict_path, "rb") as handle:
            attr_dict = pickle.load(handle)
        if "init_params_" not in attr_dict.keys():
            raise ValueError(
                "No init_params_ were saved by the model. Check out the developers guide if creating custom models."
            )
        # get the parameters for the class init signiture
        init_params = attr_dict.pop("init_params_")
        # grab all the parameters execept for kwargs (is a dict)
        non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
        # expand out kwargs
        kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        model = cls(adata, **non_kwargs, **kwargs)
        for attr, val in attr_dict.items():
            setattr(model, attr, val)
        use_cuda = use_cuda and torch.cuda.is_available()

        if use_cuda:
            model.model.load_state_dict(torch.load(model_path))
            model.model.cuda()
        else:
            device = torch.device("cpu")
            model.model.load_state_dict(torch.load(model_path, map_location=device))
        model.model.eval()
        model._validate_anndata(adata)
        return model

    def __repr__(
        self,
    ):
        summary_string = self._model_summary_string + "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        return summary_string
