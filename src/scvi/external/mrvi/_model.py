from __future__ import annotations

import os
import warnings
from typing import TYPE_CHECKING, Literal

import torch

from scvi import settings
from scvi.data._constants import _MODEL_NAME_KEY
from scvi.external.mrvi_torch import TorchMRVI
from scvi.model.base import BaseMinifiedModeModelClass
from scvi.model.base._constants import SAVE_KEYS
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from anndata import AnnData
    from lightning import LightningDataModule

    from scvi._types import AnnOrMuData


Backend = Literal["torch", "jax"]


class MRVI(BaseMinifiedModeModelClass):
    """
    Multi-resolution Variational Inference (MrVI).

    This is a convenience wrapper that instantiates the Torch or JAX
    implementation based on `backend` and returns that instance.

    Parameters
    ----------
    adata
        AnnData object that has been registered via the appropriate `setup_anndata`.
    backend
        Which backend to use: "torch" or "jax".
    registry
        (Torch-only) Registry dict for loading from saved state.
    **model_kwargs
        Extra keyword args forwarded to the selected implementation.

    Notes
    -----
    - When `backend="torch"`, this returns an instance of `TorchMRVI`.
    - When `backend="jax"`, this returns an instance of `JaxMRVI`.
    """

    def __new__(
        cls,
        adata=None,
        *,
        backend: Backend = "jax",
        registry: dict | None = None,
        **model_kwargs,
    ):
        backend = backend.lower()
        cls.backend = backend
        if backend == "torch":
            return TorchMRVI(adata=adata, registry=registry, **model_kwargs)
        elif backend == "jax":
            from scvi.external.mrvi_jax import JaxMRVI

            return JaxMRVI(adata=adata, **model_kwargs)
        else:
            raise ValueError(f"Unknown backend '{backend}'. Use 'torch' or 'jax'.")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        sample_key: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        backend: Backend = "jax",
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_sample_key)s
        %(param_batch_key)s
        %(param_labels_key)s
        backend
            Which backend to use: "torch" or "jax".
        **kwargs
            Additional keyword arguments passed into
            :meth:`~scvi.data.AnnDataManager.register_fields`.
        """
        backend = backend.lower()
        if backend == "torch":
            warnings.warn(
                "MRVI model is being setup with PyTorch backend",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            TorchMRVI.setup_anndata(
                adata=adata,
                layer=layer,
                sample_key=sample_key,
                batch_key=batch_key,
                labels_key=labels_key,
                **kwargs,
            )
        elif backend == "jax":
            from scvi.external.mrvi_jax import JaxMRVI

            warnings.warn(
                "MRVI model is being setup with JAX backend",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            JaxMRVI.setup_anndata(
                adata=adata,
                layer=layer,
                sample_key=sample_key,
                batch_key=batch_key,
                labels_key=labels_key,
                **kwargs,
            )
        else:
            raise ValueError(f"Unknown backend '{backend}'. Use 'torch' or 'jax'.")

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
        if cls.backend == "torch":
            warnings.warn(
                "MRVI model is being loaded with PyTorch backend",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

            registry = peek_loaded_model_registry(dir_path, prefix)
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] == "JaxMRVI":
                raise ValueError(
                    "It appears you are trying to load a TORCH MRVI model with a JAX MRVI model"
                )
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] == "MRVI":
                raise ValueError(
                    "It appears you are trying to load a TORCH MRVI model "
                    "with a previous version JAX MRVI model"
                )

            return TorchMRVI.load(
                dir_path,
                adata=adata,
                accelerator=accelerator,
                device=device,
                prefix=prefix,
                backup_url=backup_url,
                datamodule=datamodule,
            )
        elif cls.backend == "jax":
            from scvi.external.mrvi_jax import JaxMRVI

            warnings.warn(
                "MRVI model is being loaded with JAX backend",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

            registry = peek_loaded_model_registry(dir_path, prefix)
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] == "TorchMRVI":
                raise ValueError(
                    "It appears you are trying to load a JAX MRVI model with a Torch MRVI model"
                )

            return JaxMRVI.load(
                dir_path,
                adata=adata,
                accelerator=accelerator,
                device=device,
                prefix=prefix,
                backup_url=backup_url,
                datamodule=datamodule,
                allowed_classes_names_list=[
                    "MRVI"
                ],  # allowing old JAX MRVI models to be loaded TODO: need to change in v1.5
            )
        else:
            raise ValueError(f"Unknown backend '{cls.backend}'. Use 'torch' or 'jax'.")


def peek_loaded_model_registry(dir_path, prefix):
    """Getting the loaded model registry to give better warnings for loading MRVI"""
    file_name_prefix = prefix or ""
    model_file_name = f"{file_name_prefix}{SAVE_KEYS.MODEL_FNAME}"
    model_path = os.path.join(dir_path, model_file_name)
    attr_dict = torch.load(model_path, weights_only=False).get(SAVE_KEYS.ATTR_DICT_KEY)
    registry = attr_dict.pop("registry_")
    return registry
