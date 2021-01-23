import logging
from typing import Union

import torch
from anndata import AnnData

from scvi.compose import FCLayers
from scvi.data import transfer_anndata_setup

from ._base_model import BaseModelClass
from ._utils import _initialize_model, _load_saved_files, _validate_var_names

logger = logging.getLogger(__name__)


class ArchesMixin:
    @classmethod
    def load_query_data(
        cls,
        adata: AnnData,
        reference_model: Union[str, BaseModelClass],
        inplace_subset_query_vars: bool = False,
        use_gpu: bool = True,
        unfrozen: bool = False,
        freeze_dropout: bool = False,
        freeze_expression: bool = True,
        freeze_decoder_first_layer: bool = True,
        freeze_batchnorm_encoder: bool = True,
        freeze_batchnorm_decoder: bool = False,
        freeze_classifier: bool = True,
    ):
        """
        Online update of a reference model with scArches algorithm [Lotfollahi20]_.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run :func:`~scvi.data.setup_anndata`,
            as AnnData is validated against the saved `scvi` setup dictionary.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        inplace_subset_query_vars
            Whether to subset and rearrange query vars inplace based on vars used to
            train reference model.
        use_gpu
            Whether to load model on GPU.
        unfrozen
            Override all other freeze options for a fully unfrozen model
        freeze_dropout
            Whether to freeze dropout during training
        freeze_expression
            Freeze neurons corersponding to expression in first layer
        freeze_decoder_first_layer
            Freeze neurons corersponding to first layer in decoder
        freeze_batchnorm_encoder
            Whether to freeze batchnorm weight and bias during training for encoder
        freeze_batchnorm_decoder
            Whether to freeze batchnorm weight and bias during training for decoder
        freeze_classifier
            Whether to freeze classifier completely. Only applies to `SCANVI`.
        """
        use_gpu = use_gpu and torch.cuda.is_available()

        if isinstance(reference_model, str):
            map_location = torch.device("cuda") if use_gpu is True else None
            (
                scvi_setup_dict,
                attr_dict,
                var_names,
                load_state_dict,
                _,
            ) = _load_saved_files(
                reference_model, load_adata=False, map_location=map_location
            )
        else:
            attr_dict = reference_model._get_user_attributes()
            attr_dict = {a[0]: a[1] for a in attr_dict if a[0][-1] == "_"}
            scvi_setup_dict = attr_dict.pop("scvi_setup_dict_")
            var_names = reference_model.adata.var_names
            load_state_dict = reference_model.module.state_dict().copy()

        if inplace_subset_query_vars:
            logger.debug("Subsetting query vars to reference vars.")
            adata._inplace_subset_var(var_names)
        _validate_var_names(adata, var_names)

        if scvi_setup_dict["scvi_version"] < "0.8":
            logger.warning(
                "Query integration should be performed using models trained with version >= 0.8"
            )

        transfer_anndata_setup(scvi_setup_dict, adata, extend_categories=True)

        model = _initialize_model(cls, adata, attr_dict, use_gpu)

        # set saved attrs for loaded model
        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        if use_gpu:
            model.module.cuda()

        # model tweaking
        new_state_dict = model.module.state_dict()
        for key, load_ten in load_state_dict.items():
            new_ten = new_state_dict[key]
            if new_ten.size() == load_ten.size():
                continue
            # new categoricals changed size
            else:
                dim_diff = new_ten.size()[-1] - load_ten.size()[-1]
                fixed_ten = torch.cat([load_ten, new_ten[..., -dim_diff:]], dim=-1)
                load_state_dict[key] = fixed_ten

        model.module.load_state_dict(load_state_dict)
        model.module.eval()

        _set_params_online_update(
            model.module,
            unfrozen=unfrozen,
            freeze_decoder_first_layer=freeze_decoder_first_layer,
            freeze_batchnorm_encoder=freeze_batchnorm_encoder,
            freeze_batchnorm_decoder=freeze_batchnorm_decoder,
            freeze_dropout=freeze_dropout,
            freeze_expression=freeze_expression,
            freeze_classifier=freeze_classifier,
        )
        model.is_trained_ = False

        return model


def _set_params_online_update(
    module,
    unfrozen,
    freeze_decoder_first_layer,
    freeze_batchnorm_encoder,
    freeze_batchnorm_decoder,
    freeze_dropout,
    freeze_expression,
    freeze_classifier,
):
    """Freeze parts of network for scArches."""
    # do nothing if unfrozen
    if unfrozen:
        return

    mod_no_grad = set(["encoder_z2_z1", "decoder_z1_z2"])
    mod_no_hooks_yes_grad = set(["l_encoder"])
    if not freeze_classifier:
        mod_no_hooks_yes_grad.add("classifier")
    parameters_yes_grad = set(["background_pro_alpha", "background_pro_log_beta"])

    def no_hook_cond(key):
        one = (not freeze_expression) and "encoder" in key
        two = (not freeze_decoder_first_layer) and "px_decoder" in key
        return one or two

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

    for key, mod in module.named_modules():
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

    for key, par in module.named_parameters():
        if requires_grad(key):
            par.requires_grad = True
        else:
            par.requires_grad = False
