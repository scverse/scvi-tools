import logging
import os
import pickle

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
        freeze_dropout: bool = False,
        freeze_batchnorm: bool = False,
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
        freeze_dropout
            Whether to freeze dropout during training
        freeze_batchnorm
            Whether to freeze batchnorm statistics during training
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
        setattr(model, "use_cuda", use_cuda)
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

        _set_params_online_update(
            model.model,
            freeze_batchnorm=freeze_batchnorm,
            freeze_dropout=freeze_dropout,
        )
        model.is_trained_ = False

        return model


def _set_params_online_update(model, freeze_batchnorm, freeze_dropout):

    mod_no_hooks = ["classifier", "encoder_z2_z1", "decoder_z1_z2"]

    for key, mod in model.named_modules():
        if isinstance(mod, FCLayers):
            mod.set_online_update_hooks()
        if isinstance(mod, torch.nn.Dropout):
            if freeze_dropout:
                mod.p = 0
        if isinstance(mod, torch.nn.BatchNorm1d):
            if freeze_batchnorm:
                mod.affine = False
                mod.track_running_stats = False

    for key, par in model.named_parameters():
        # gets the linear layer
        if "fc_layers" in key and ".0." in key:
            for m in mod_no_hooks:
                if m in key:
                    par.requires_grad = False
                else:
                    par.requires_grad = True
        else:
            par.requires_grad = False
