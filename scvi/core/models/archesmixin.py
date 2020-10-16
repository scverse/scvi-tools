import logging

import torch
from anndata import AnnData

from scvi.core.modules._base import FCLayers
from scvi.data import transfer_anndata_setup
from scvi.core.models._utils import (
    _initialize_model,
    _load_saved_files,
    _validate_var_names,
)

logger = logging.getLogger(__name__)


class ArchesMixin:
    @classmethod
    def load_query_data(
        cls,
        adata: AnnData,
        dir_path: str,
        use_cuda: bool = True,
        freeze_dropout: bool = False,
        freeze_expression: bool = True,
        freeze_batchnorm_encoder: bool = False,
        freeze_batchnorm_decoder: bool = True,
    ):
        """
        Online update of a reference model with scArches algorithm [Lotfollahi20]_.

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
        freeze_expression
            Freeze neurons corersponding to expression in first layer
        freeze_batchnorm_encoder
            Whether to freeze batchnorm weight and bias during training for encoder
        freeze_batchnorm_decoder
            Whether to freeze batchnorm weight and bias during training for decoder
        """
        use_cuda = use_cuda and torch.cuda.is_available()
        map_location = torch.device("cpu") if use_cuda is False else None
        scvi_setup_dict, attr_dict, var_names, load_state_dict, _ = _load_saved_files(
            dir_path, load_adata=False, map_location=map_location
        )

        _validate_var_names(adata, var_names)
        transfer_anndata_setup(scvi_setup_dict, adata, extend_categories=True)
        model = _initialize_model(cls, adata, attr_dict, use_cuda)

        # set saved attrs for loaded model
        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        if use_cuda:
            model.model.cuda()

        # model tweaking
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
            freeze_batchnorm_encoder=freeze_batchnorm_encoder,
            freeze_batchnorm_decoder=freeze_batchnorm_decoder,
            freeze_dropout=freeze_dropout,
            freeze_expression=freeze_expression,
        )
        model.is_trained_ = False

        return model


def _set_params_online_update(
    model,
    freeze_batchnorm_encoder,
    freeze_batchnorm_decoder,
    freeze_dropout,
    freeze_expression,
):
    """Freeze parts of network for scArches."""
    mod_no_grad = set(["encoder_z2_z1", "decoder_z1_z2"])
    mod_no_hooks_yes_grad = set(["l_encoder"])
    parameters_yes_grad = set(["background_pro_alpha", "background_pro_log_beta"])

    def no_hook_cond(key):
        return (not freeze_expression) and "encoder" in key

    def requires_grad(key):
        mod_name = key.split(".")[0]
        # linear weights and bias that need grad
        one = "fc_layers" in key and ".0." in key and mod_name not in mod_no_grad
        # modules that need grad
        two = mod_name in mod_no_hooks_yes_grad
        three = sum([p in key for p in parameters_yes_grad]) > 0
        # batch norm option
        four = (
            "fc_layers" in key
            and ".1." in key
            and "encoder" in key
            and (not freeze_batchnorm_encoder)
        )
        five = (
            "fc_layers" in key
            and ".1." in key
            and "decoder" in key
            and (not freeze_batchnorm_decoder)
        )
        if one or two or three or four or five:
            return True
        else:
            return False

    for key, mod in model.named_modules():
        # skip over protected modules
        if key.split(".")[0] in mod_no_hooks_yes_grad:
            continue
        if isinstance(mod, FCLayers):
            hook_first_layer = False if no_hook_cond(key) else True
            mod.set_online_update_hooks(hook_first_layer)
        if isinstance(mod, torch.nn.Dropout):
            if freeze_dropout:
                mod.p = 0
        # momentum freezes the running stats of batchnorm
        freeze_batchnorm = ("decoder" in key and freeze_batchnorm_decoder) or (
            "encoder" in key and freeze_batchnorm_encoder
        )
        if isinstance(mod, torch.nn.BatchNorm1d) and freeze_batchnorm:
            mod.momentum = 0

    for key, par in model.named_parameters():
        if requires_grad(key):
            par.requires_grad = True
        else:
            par.requires_grad = False
