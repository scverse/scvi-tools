import logging
import os
import pickle
from typing import Optional

import torch
from anndata import AnnData

from scvi.core.modules._base import FCLayers
from scvi.data import transfer_anndata_setup

logger = logging.getLogger(__name__)


class ArchesMixin:
    @classmethod
    def load_query_data(
        cls,
        adata: AnnData,
        dir_path: str,
        use_cuda: bool = True,
        n_epochs: Optional[int] = None,
        train_size: float = 0.9,
        test_size: Optional[float] = None,
        lr: float = 1e-3,
        n_epochs_kl_warmup: int = 400,
        n_iter_kl_warmup: Optional[int] = None,
        frequency: Optional[int] = None,
        train_fun_kwargs: dict = {},
        **kwargs,
    ):
        """
        Online update of a reference model using method of scArches.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :func:`~scvi.data.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
        dir_path
            Path to saved outputs for reference model.
        use_cuda
            Whether to load model on GPU.
        n_epochs
            Number of passes through the query dataset.
        train_size
            Size of training set in the range [0.0, 1.0].
        test_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + test_size < 1`, the remaining cells belong to a validation set.
        lr
            Learning rate for optimization.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        n_iter_kl_warmup
            Number of minibatches for scaling term on KL divergence to go from 0 to 1.
            To use, set to not `None` and set `n_epochs_kl_warmup` to `None`.
        frequency
            Frequency with which metrics are computed on the data for train/test/val sets.
        train_fun_kwargs
            Keyword args for the train method of :class:`~scvi.core.trainers.UnsupervisedTrainer`.
        **kwargs
            Other keyword args for :class:`~scvi.core.trainers.UnsupervisedTrainer`.
        """
        model_path = os.path.join(dir_path, "model_params.pt")
        setup_dict_path = os.path.join(dir_path, "attr.pkl")
        with open(setup_dict_path, "rb") as handle:
            attr_dict = pickle.load(handle)

        transfer_anndata_setup(
            attr_dict["scvi_setup_dict_"], adata, extend_categories=True
        )
        attr_dict["scvi_setup_dict_"] = adata.uns["_scvi"]

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
            model.model.cuda()

        # model tweaking
        load_state_dict = torch.load(
            model_path, map_location=torch.device("cpu") if not use_cuda else None
        )
        new_state_dict = model.model.state_dict()
        for key, load_ten in load_state_dict.items():
            new_ten = new_state_dict[key]
            if new_ten.size() == load_ten.size():
                continue
            # new categoricals changed size
            else:
                dim_diff = new_ten.size()[-1] - load_ten.size()[-1]
                fixed_ten = torch.cat([load_ten, new_ten[..., -dim_diff:]], dim=-1)
                load_state_dict[key] = fixed_ten

        model.model.load_state_dict(load_state_dict)
        model.model.eval()

        _set_params_online_update(model.model)
        model.is_trained_ = False

        return model


def _set_params_online_update(model):

    for key, mod in model.named_modules():
        if isinstance(mod, FCLayers):
            mod.set_online_update_hooks()
        if isinstance(mod, torch.nn.Dropout):
            mod.p = 0
        if isinstance(mod, torch.nn.BatchNorm1d):
            mod.affine = False
            mod.track_running_stats = False

    for key, par in model.named_parameters():
        # gets the linear layer
        if "fc_layers" in key and ".0." in key:
            par.requires_grad = True
        else:
            par.requires_grad = False
