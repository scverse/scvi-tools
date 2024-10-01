import logging
import warnings
from copy import deepcopy

import anndata
import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from mudata import MuData
from scipy.sparse import csr_matrix

from scvi import settings
from scvi._types import AnnOrMuData
from scvi.data import _constants
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY, _SETUP_METHOD_NAME
from scvi.model._utils import parse_device_args
from scvi.model.base._save_load import (
    _get_var_names,
    _initialize_model,
    _load_saved_files,
    _validate_var_names,
)
from scvi.nn import FCLayers
from scvi.utils._docstrings import devices_dsp

from ._base_model import BaseModelClass

logger = logging.getLogger(__name__)

MIN_VAR_NAME_RATIO = 0.8


class ArchesMixin:
    """Universal scArches implementation."""

    @classmethod
    @devices_dsp.dedent
    def load_query_data(
        cls,
        adata: AnnOrMuData,
        reference_model: str | BaseModelClass,
        inplace_subset_query_vars: bool = False,
        accelerator: str = "auto",
        device: int | str = "auto",
        unfrozen: bool = False,
        freeze_dropout: bool = False,
        freeze_expression: bool = True,
        freeze_decoder_first_layer: bool = True,
        freeze_batchnorm_encoder: bool = True,
        freeze_batchnorm_decoder: bool = False,
        freeze_classifier: bool = True,
    ):
        """Online update of a reference model with scArches algorithm :cite:p:`Lotfollahi21`.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        inplace_subset_query_vars
            Whether to subset and rearrange query vars inplace based on vars used to
            train reference model.
        %(param_accelerator)s
        %(param_device)s
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
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        attr_dict, var_names, load_state_dict = _get_loaded_data(reference_model, device=device)

        if isinstance(adata, MuData):
            for modality in adata.mod:
                if inplace_subset_query_vars:
                    logger.debug(f"Subsetting {modality} query vars to reference vars.")
                    adata[modality]._inplace_subset_var(var_names[modality])
                _validate_var_names(adata[modality], var_names[modality])

        else:
            if inplace_subset_query_vars:
                logger.debug("Subsetting query vars to reference vars.")
                adata._inplace_subset_var(var_names)
            _validate_var_names(adata, var_names)

        if inplace_subset_query_vars:
            logger.debug("Subsetting query vars to reference vars.")
            adata._inplace_subset_var(var_names)
        _validate_var_names(adata, var_names)

        registry = attr_dict.pop("registry_")
        if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
            raise ValueError("It appears you are loading a model from a different class.")

        if _SETUP_ARGS_KEY not in registry:
            raise ValueError(
                "Saved model does not contain original setup inputs. "
                "Cannot load the original setup."
            )

        setup_method = getattr(cls, registry[_SETUP_METHOD_NAME])
        setup_method(
            adata,
            source_registry=registry,
            extend_categories=True,
            allow_missing_labels=True,
            **registry[_SETUP_ARGS_KEY],
        )

        model = _initialize_model(cls, adata, attr_dict)
        adata_manager = model.get_anndata_manager(adata, required=True)

        version_split = adata_manager.registry[_constants._SCVI_VERSION_KEY].split(".")
        if int(version_split[1]) < 8 and int(version_split[0]) == 0:
            warnings.warn(
                "Query integration should be performed using models trained with "
                "version >= 0.8",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        model.to_device(device)

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

    @staticmethod
    def prepare_query_anndata(
        adata: AnnData,
        reference_model: str | BaseModelClass,
        return_reference_var_names: bool = False,
        inplace: bool = True,
    ) -> AnnData | pd.Index | None:
        """Prepare data for query integration.

        This function will return a new AnnData object with padded zeros
        for missing features, as well as correctly sorted features.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        return_reference_var_names
            Only load and return reference var names if True.
        inplace
            Whether to subset and rearrange query vars inplace or return new AnnData.

        Returns
        -------
        Query adata ready to use in `load_query_data` unless `return_reference_var_names`
        in which case a pd.Index of reference var names is returned.
        """
        _, var_names, _ = _get_loaded_data(reference_model, device="cpu")
        var_names = pd.Index(var_names)

        if return_reference_var_names:
            return var_names

        return _pad_and_sort_query_anndata(adata, var_names, inplace)

    @staticmethod
    def prepare_query_mudata(
        mdata: MuData,
        reference_model: str | BaseModelClass,
        return_reference_var_names: bool = False,
        inplace: bool = True,
    ) -> MuData | dict[str, pd.Index] | None:
        """Prepare multimodal dataset for query integration.

        This function will return a new MuData object such that the
        AnnData objects for individual modalities are given padded zeros
        for missing features, as well as correctly sorted features.

        Parameters
        ----------
        mdata
            MuData organized in the same way as data used to train model.
            It is not necessary to run setup_mudata,
            as MuData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        return_reference_var_names
            Only load and return reference var names if True.
        inplace
            Whether to subset and rearrange query vars inplace or return new MuData.

        Returns
        -------
        Query mudata ready to use in `load_query_data` unless `return_reference_var_names`
        in which case a dictionary of pd.Index of reference var names is returned.
        """
        attr_dict, var_names, _ = _get_loaded_data(reference_model, device="cpu")

        for modality in var_names.keys():
            var_names[modality] = pd.Index(var_names[modality])

        if return_reference_var_names:
            return var_names

        reference_modalities_dict = attr_dict["registry_"][_SETUP_ARGS_KEY]["modalities"]

        reference_modalities = reference_modalities_dict.values()
        query_modalities = mdata.mod

        for mod in reference_modalities:
            if mod not in query_modalities:
                raise ValueError(
                    "Query MuData does not contain same modalities as reference. "
                    "Cannot load the original setup."
                )

        adata_dict = {}
        for modality in reference_modalities:
            adata_out = _pad_and_sort_query_anndata(
                adata=mdata[modality],
                reference_var_names=var_names[modality],
                inplace=inplace,
            )
            adata_dict[modality] = adata_out

        if not inplace:
            return MuData(adata_dict)


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

    mod_inference_mode = {"encoder_z2_z1", "decoder_z1_z2"}
    mod_no_hooks_yes_grad = {"l_encoder"}
    if not freeze_classifier:
        mod_no_hooks_yes_grad.add("classifier")
    parameters_yes_grad = {"background_pro_alpha", "background_pro_log_beta"}

    def no_hook_cond(key):
        one = (not freeze_expression) and "encoder" in key
        two = (not freeze_decoder_first_layer) and "px_decoder" in key
        return one or two

    def requires_grad(key):
        mod_name = key.split(".")[0]
        # linear weights and bias that need grad
        one = "fc_layers" in key and ".0." in key and mod_name not in mod_inference_mode
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


def _get_loaded_data(reference_model, device=None):
    if isinstance(reference_model, str):
        attr_dict, var_names, load_state_dict, _ = _load_saved_files(
            reference_model, load_adata=False, map_location=device
        )
    else:
        attr_dict = reference_model._get_user_attributes()
        attr_dict = {a[0]: a[1] for a in attr_dict if a[0][-1] == "_"}
        var_names = _get_var_names(reference_model.adata)
        load_state_dict = deepcopy(reference_model.module.state_dict())

    return attr_dict, var_names, load_state_dict


def _pad_and_sort_query_anndata(
    adata: AnnData,
    reference_var_names: pd.Index,
    inplace: bool,
):
    intersection = adata.var_names.intersection(reference_var_names)
    inter_len = len(intersection)
    if inter_len == 0:
        raise ValueError(
            "No reference var names found in query data. "
            "Please rerun with return_reference_var_names=True "
            "to see reference var names."
        )

    ratio = inter_len / len(reference_var_names)
    logger.info(f"Found {ratio * 100}% reference vars in query data.")
    if ratio < MIN_VAR_NAME_RATIO:
        warnings.warn(
            f"Query data contains less than {MIN_VAR_NAME_RATIO:.0%} of reference "
            "var names. This may result in poor performance.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
    genes_to_add = reference_var_names.difference(adata.var_names)
    needs_padding = len(genes_to_add) > 0
    if needs_padding:
        padding_mtx = csr_matrix(np.zeros((adata.n_obs, len(genes_to_add))))
        adata_padding = AnnData(
            X=padding_mtx.copy(),
            layers={layer: padding_mtx.copy() for layer in adata.layers},
        )
        adata_padding.var_names = genes_to_add
        adata_padding.obs_names = adata.obs_names
        # Concatenate object
        adata_out = anndata.concat(
            [adata, adata_padding],
            axis=1,
            join="outer",
            index_unique=None,
            merge="unique",
        )
    else:
        adata_out = adata

    # also covers the case when new adata has more var names than old
    if not reference_var_names.equals(adata_out.var_names):
        adata_out._inplace_subset_var(reference_var_names)

    if inplace:
        if adata_out is not adata:
            adata._init_as_actual(adata_out)
    else:
        return adata_out
