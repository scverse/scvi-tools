from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import scvi
import torch
from lightning import LightningDataModule
from scvi import REGISTRY_KEYS
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY, _SETUP_METHOD_NAME
from scvi.model._utils import parse_device_args
from scvi.model.base import BaseModelClass
from scvi.model.base._archesmixin import ArchesMixin, _get_loaded_data, _initialize_model, _validate_var_names
from torch import nn

from scvi.external.drvi.nn import FCLayers

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

logger = logging.getLogger(__name__)


class DRVIArchesMixin(ArchesMixin):
    """Universal scArches implementation."""

    @classmethod
    def load_query_data(
        cls,
        adata: AnnData = None,
        reference_model: str | BaseModelClass = None,
        registry: dict = None,
        inplace_subset_query_vars: bool = False,
        accelerator: str = "auto",
        device: int | str = "auto",
        unfrozen: bool = False,
        freeze_dropout: bool = False,
        freeze_shared_emb: bool = True,
        freeze_encoder: bool = True,
        freeze_decoder: bool = True,
        reset_encoder: bool = False,
        reset_decoder: bool = False,
        freeze_batchnorm_encoder: bool = True,
        freeze_batchnorm_decoder: bool = False,
        datamodule: LightningDataModule | None = None,
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
        freeze_shared_emb
            Whether to freeze shared embeddings if any
        freeze_encoder
            Whether to freeze encoder
        freeze_decoder
            Whether to freeze decoder
        freeze_batchnorm_encoder
            Whether to freeze encoder batchnorms' weight and bias during transfer
        freeze_batchnorm_decoder
            Whether to freeze decoder batchnorms' weight and bias during transfer
        """
        if reference_model is None:
            raise ValueError("Please provide a reference model as string or loaded model.")
        if adata is None and registry is None:
            raise ValueError("Please provide either an AnnData or a registry dictionary.")

        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        # We limit to [:3] as from scvi version 1.1.5 additional output (pyro_param_store) is returned
        attr_dict, var_names, load_state_dict = _get_loaded_data(reference_model, device=device)[:3]

        if adata:
            if inplace_subset_query_vars:
                logger.debug("Subsetting query vars to reference vars.")
                adata._inplace_subset_var(var_names)
            _validate_var_names(adata, var_names)

            registry = attr_dict.pop("registry_")
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
                raise ValueError("It appears you are loading a model from a different class.")

            if _SETUP_ARGS_KEY not in registry:
                raise ValueError("Saved model does not contain original setup inputs. Cannot load the original setup.")

            if registry[_SETUP_METHOD_NAME] != "setup_datamodule":
                setup_method = getattr(cls, registry[_SETUP_METHOD_NAME])
                setup_method(
                    adata,
                    source_registry=registry,
                    extend_categories=True,
                    allow_missing_labels=True,
                    **registry[_SETUP_ARGS_KEY],
                )

            cls.setup_anndata(
                adata,
                source_registry=registry,
                extend_categories=True,
                allow_missing_labels=True,
                **registry[_SETUP_ARGS_KEY],
            )

        model = _initialize_model(cls, adata, registry, attr_dict, datamodule)
        adata_manager = model.get_anndata_manager(adata, required=True)

        previous_n_batch = registry["field_registries"][REGISTRY_KEYS.BATCH_KEY]["summary_stats"]["n_batch"]
        n_batch = model.summary_stats.n_batch
        if REGISTRY_KEYS.CAT_COVS_KEY in adata_manager.data_registry:
            previous_n_cats_per_cov = [previous_n_batch] + list(
                registry["field_registries"][REGISTRY_KEYS.CAT_COVS_KEY]["state_registry"]["n_cats_per_key"]
            )
            n_cats_per_cov = [n_batch] + list(
                model.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            )
        else:
            previous_n_cats_per_cov = [previous_n_batch]
            n_cats_per_cov = [n_batch]

        model.to_device(device)

        reloaded_tensor_keys = []
        # model tweaking
        new_state_dict = model.module.state_dict()
        for key, load_ten in load_state_dict.items():
            new_ten = new_state_dict[key]
            if reset_decoder and "decoder" in key:
                assert not freeze_decoder
                print(f"Resetting {key} in decoder")
                load_state_dict[key] = new_ten
                reloaded_tensor_keys.append(key)
                continue
            if reset_encoder and "z_encoder" in key:
                assert not freeze_encoder
                print(f"Resetting {key} in encoder")
                load_state_dict[key] = new_ten
                reloaded_tensor_keys.append(key)
                continue
            if new_ten.size() == load_ten.size():
                continue
            # new categoricals changed size
            else:
                print(f"Resizing {key} from {load_ten.size()} to {new_ten.size()}")
                if "emb_list" in key:
                    # Extend embeddings along dim 0 (n_emb)
                    dim_diff = new_ten.size()[0] - load_ten.size()[0]
                    assert new_ten.size()[1] == load_ten.size()[1]
                    fixed_ten = torch.cat([load_ten, new_ten[-dim_diff:, :]], dim=0)
                    load_state_dict[key] = fixed_ten
                    reloaded_tensor_keys.append(key)
                else:
                    if new_ten.dim() == load_ten.dim() == 2:
                        # 2D tensors
                        # Extend linear layers along dim 1 (input)
                        assert (
                            sum(n_cats_per_cov) - sum(previous_n_cats_per_cov)
                            == new_ten.size()[-1] - load_ten.size()[-1]
                        )
                        assert new_ten.size()[0] == load_ten.size()[0]
                        # Keep data flow (normal nodes) as is
                        cum_n_cat_old = new_ten.size()[1] - sum(n_cats_per_cov)
                        cum_n_cat_new = cum_n_cat_old
                        fixed_ten = [load_ten[:, :cum_n_cat_old]] if cum_n_cat_old > 0 else []
                        # Iterate and get old covariates deom load_ten and new ones from new_ten
                        for n_cat_new, n_cat_old in zip(n_cats_per_cov, previous_n_cats_per_cov, strict=False):
                            fixed_ten.append(load_ten[:, cum_n_cat_old : cum_n_cat_old + n_cat_old])
                            if n_cat_new > n_cat_old:
                                fixed_ten.append(new_ten[:, cum_n_cat_new : cum_n_cat_new + n_cat_new - n_cat_old])
                            cum_n_cat_old += n_cat_old
                            cum_n_cat_new += n_cat_new
                        # Concat and set as init tensor
                        fixed_ten = torch.cat([t.to(device) for t in fixed_ten], dim=1)
                        load_state_dict[key] = fixed_ten
                        reloaded_tensor_keys.append(key)
                    elif new_ten.dim() == load_ten.dim() == 3:
                        # 3D tensors
                        # extend weight of stacked linears along dim 1 (input)
                        assert (
                            sum(n_cats_per_cov) - sum(previous_n_cats_per_cov) == new_ten.size()[1] - load_ten.size()[1]
                        )
                        assert new_ten.size()[0] == load_ten.size()[0]
                        assert new_ten.size()[2] == load_ten.size()[2]
                        # Keep data flow (normal nodes) as is
                        cum_n_cat_old = new_ten.size()[1] - sum(n_cats_per_cov)
                        cum_n_cat_new = cum_n_cat_old
                        fixed_ten = [load_ten[:, :cum_n_cat_old, :]] if cum_n_cat_old > 0 else []
                        # Iterate and get old covariates deom load_ten and new ones from new_ten
                        for n_cat_new, n_cat_old in zip(n_cats_per_cov, previous_n_cats_per_cov, strict=False):
                            fixed_ten.append(load_ten[:, cum_n_cat_old : cum_n_cat_old + n_cat_old, :])
                            if n_cat_new > n_cat_old:
                                fixed_ten.append(new_ten[:, cum_n_cat_new : cum_n_cat_new + n_cat_new - n_cat_old])
                            cum_n_cat_old += n_cat_old
                            cum_n_cat_new += n_cat_new
                        # Concat and set as init tensor
                        fixed_ten = torch.cat([t.to(device) for t in fixed_ten], dim=1)
                        load_state_dict[key] = fixed_ten
                        reloaded_tensor_keys.append(key)
                    else:
                        raise NotImplementedError()

        model.module.load_state_dict(load_state_dict)
        model.module.eval()

        # Make hooks to make gradients zero while training
        _set_params_online_update(
            model.module,
            reloaded_tensor_keys,
            unfrozen=unfrozen,
            previous_n_cats_per_cov=previous_n_cats_per_cov,
            n_cats_per_cov=n_cats_per_cov,
            freeze_dropout=freeze_dropout,
            freeze_shared_emb=freeze_shared_emb,
            freeze_encoder=freeze_encoder,
            freeze_decoder=freeze_decoder,
            freeze_batchnorm_encoder=freeze_batchnorm_encoder,
            freeze_batchnorm_decoder=freeze_batchnorm_decoder,
        )
        model.is_trained_ = False

        return model


def _set_params_online_update(
    module: nn.Module,
    reloaded_tensor_keys: list[str],
    unfrozen: bool,
    previous_n_cats_per_cov: Sequence[int] | None,
    n_cats_per_cov: Sequence[int] | None,
    freeze_dropout: bool,
    freeze_shared_emb: bool,
    freeze_encoder: bool,
    freeze_decoder: bool,
    freeze_batchnorm_encoder: bool,
    freeze_batchnorm_decoder: bool,
) -> None:
    """Freeze parts of network for scArches."""
    # do nothing if unfrozen
    if unfrozen:
        return

    if freeze_shared_emb:
        if hasattr(module, "shared_covariate_emb") and module.shared_covariate_emb is not None:
            print(f"Freezing top {previous_n_cats_per_cov} items in Shared Emb.")
            assert tuple(module.shared_covariate_emb.num_embeddings) == tuple(n_cats_per_cov)
            module.shared_covariate_emb.freeze_top_embs(previous_n_cats_per_cov)

    def requires_grad(key: str) -> bool:
        if not freeze_decoder and "decoder" in key:
            return True
        if not freeze_encoder and "z_encoder" in key:
            return True
        # Extended tensors should have gradients but we put some gradients to zero by online hooks
        if key in reloaded_tensor_keys:
            return True
        return False

    # General changes
    for key, par in module.named_parameters():
        if requires_grad(key):
            par.requires_grad = True
        else:
            print(f"Freezing key {key}")
            par.requires_grad = False

    # specific changes
    for key, mod in module.named_modules():
        # skip over protected modules
        if isinstance(mod, FCLayers):
            freeze_fc_layers = ("decoder" in key and freeze_decoder) or ("z_encoder" in key and freeze_encoder)
            # This will make requires_frad for layers after one_hot to True
            if freeze_fc_layers:
                print(f"Setting hooks for {key}")
                mod.set_online_update_hooks(previous_n_cats_per_cov, n_cats_per_cov)
        elif isinstance(mod, nn.Dropout):
            if freeze_dropout:
                print(f"Freezing dropout {key}")
                mod.p = 0
        elif isinstance(mod, nn.BatchNorm1d | nn.LayerNorm):
            assert hasattr(mod, "freeze")
            freeze_batchnorm = ("decoder" in key and freeze_batchnorm_decoder) or (
                "z_encoder" in key and freeze_batchnorm_encoder
            )
            if freeze_batchnorm:
                print(f"Freezing normalization layer {key}")
                mod.freeze(freeze_batchnorm)
                mod.momentum = 0.0
