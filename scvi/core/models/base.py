import inspect
import logging
import os
import pickle
from abc import ABC, abstractmethod
from typing import Optional, Sequence

import numpy as np
import rich
import torch
from anndata import AnnData
from rich.text import Text
from sklearn.model_selection._split import _validate_shuffle_split

from scvi import _CONSTANTS, settings
from scvi.core.lightning import VAETask, Trainer
from scvi.data import get_from_registry, transfer_anndata_setup
from scvi.data._utils import (
    _check_anndata_setup_equivalence,
    _check_nonnegative_integers,
)
from scvi.core.models._utils import (
    _initialize_model,
    _load_saved_files,
    _validate_var_names,
)

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
        self.history_ = None

    def _make_scvi_dl(
        self,
        adata: AnnData,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        shuffle=False,
        **data_loader_kwargs,
    ):
        """Create a ScviDataLoader object for data iteration."""
        if batch_size is None:
            batch_size = settings.batch_size
        if indices is None:
            indices = np.arange(adata.n_obs)
        post = self._scvi_dl_class(
            adata,
            shuffle=shuffle,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=batch_size,
            **data_loader_kwargs,
        )
        return post

    def _train_test_val_split(
        self, adata, train_size=0.9, validation_size=None, **kwargs
    ):
        """
        Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

        If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

        Parameters
        ----------
        train_size
            float, or None (default is 0.9)
        validation_size
            float, or None (default is None)
        """
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        n = len(adata)
        try:
            n_train, n_val = _validate_shuffle_split(n, validation_size, train_size)
        except ValueError:
            if train_size != 1.0:
                raise ValueError(
                    "Choice of train_size={} and validation_size={} not understood".format(
                        train_size, validation_size
                    )
                )
            n_train, n_val = n, 0
        random_state = np.random.RandomState(seed=settings.seed)
        permutation = random_state.permutation(n)
        indices_validation = permutation[:n_val]
        indices_train = permutation[n_val : (n_val + n_train)]
        indices_test = permutation[(n_val + n_train) :]

        return (
            self._make_scvi_dl(adata, indices=indices_train, shuffle=True, **kwargs),
            self._make_scvi_dl(
                adata, indices=indices_validation, shuffle=True, **kwargs
            ),
            self._make_scvi_dl(adata, indices=indices_test, shuffle=True, **kwargs),
        )

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

    @property
    def history(self):
        """Returns computed metrics during training."""
        return self.history_

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

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: bool = True,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        lr
            Learning rate for optimization.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        n_iter_kl_warmup
            Number of minibatches for scaling term on KL divergence to go from 0 to 1.
            To use, set to not `None` and set `n_epochs_kl_warmup` to `None`.
        **kwargs
            Other keyword args for :class:`~scvi.core.lightning.Trainer`.
        """
        if use_gpu is not None:
            use_gpu = use_gpu and torch.cuda.is_available()
        else:
            use_gpu = torch.cuda.is_available()

        if self.is_trained_ is False:
            self._pl_task = VAETask(
                self.model,
            )
            if use_gpu:
                gpus = 1
                pin_memory = True
            else:
                gpus = None
                pin_memory = False
            self.trainer = Trainer(max_epochs=max_epochs, gpus=gpus, **kwargs)
            train_dl, val_dl, test_dl = self._train_test_val_split(
                self.adata,
                train_size=train_size,
                validation_size=validation_size,
                pin_memory=pin_memory,
            )
            self.train_indices_ = train_dl.indices
            self.test_indices_ = test_dl.indices
            self.validation_indices_ = val_dl.indices

        self.trainer.fit(self._pl_task, train_dl, val_dl)
        self.model.eval()
        if use_gpu:
            self.model.cuda()
        self.is_trained_ = True

    def save(
        self,
        dir_path: str,
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
        overwrite
            Overwrite existing data or not. If `False` and directory
            already exists at `dir_path`, error will be raised.
        save_anndata
            If True, also saves the anndata
        anndata_write_kwargs
            Kwargs for anndata write function
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

        if save_anndata:
            self.adata.write(
                os.path.join(dir_path, "adata.h5ad"), **anndata_write_kwargs
            )

        model_save_path = os.path.join(dir_path, "model_params.pt")
        attr_save_path = os.path.join(dir_path, "attr.pkl")
        varnames_save_path = os.path.join(dir_path, "var_names.csv")

        var_names = self.adata.var_names.astype(str)
        var_names = var_names.to_numpy()
        np.savetxt(varnames_save_path, var_names, fmt="%s")

        torch.save(self.model.state_dict(), model_save_path)
        with open(attr_save_path, "wb") as f:
            pickle.dump(user_attributes, f)

    @classmethod
    def load(
        cls,
        dir_path: str,
        adata: Optional[AnnData] = None,
        use_cuda: bool = False,
    ):
        """
        Instantiate a model from the saved output.

        Parameters
        ----------
        dir_path
            Path to saved outputs.
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :func:`~scvi.data.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
            If None, will check for and load anndata saved with the model.
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
        load_adata = adata is None
        use_cuda = use_cuda and torch.cuda.is_available()
        map_location = torch.device("cpu") if use_cuda is False else None

        (
            scvi_setup_dict,
            attr_dict,
            var_names,
            model_state_dict,
            new_adata,
        ) = _load_saved_files(dir_path, load_adata, map_location=map_location)
        adata = new_adata if new_adata is not None else adata

        _validate_var_names(adata, var_names)
        transfer_anndata_setup(scvi_setup_dict, adata)
        model = _initialize_model(cls, adata, attr_dict, use_cuda)

        # set saved attrs for loaded model
        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        model.model.load_state_dict(model_state_dict)
        if use_cuda:
            model.model.cuda()

        model.model.eval()
        model._validate_anndata(adata)

        return model

    def __repr__(
        self,
    ):
        summary_string = self._model_summary_string
        summary_string += "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        rich.print(summary_string)

        command = "scvi.data.view_anndata_setup(model.adata)"
        command_len = len(command)
        print_adata_str = "\n\nTo print summary of associated AnnData, use: " + command
        text = Text(print_adata_str)
        text.stylize(
            "dark_violet", len(print_adata_str) - command_len, len(print_adata_str)
        )
        console = rich.console.Console()
        console.print(text)
        return ""
