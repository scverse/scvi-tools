from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

from scvi import settings
from scvi.external.mrvi_jax import JaxMRVI
from scvi.external.mrvi_torch import TorchMRVI
from scvi.model.base import BaseMinifiedModeModelClass
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
        if backend == "torch":
            return TorchMRVI(adata=adata, registry=registry, **model_kwargs)
        elif backend == "jax":
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
                DeprecationWarning,
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
            warnings.warn(
                "MRVI model is being setup with JAX backend",
                DeprecationWarning,
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
        backend: Backend = "jax",
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
        backend
            Which backend to use: "torch" or "jax".

        Returns
        -------
        Model with loaded state dictionaries.

        Examples
        --------
        >>> model = ModelClass.load(save_path, adata)
        >>> model.get_....
        """
        backend = backend.lower()
        if backend == "torch":
            warnings.warn(
                "MRVI model is being loaded with PyTorch backend",
                DeprecationWarning,
                stacklevel=settings.warnings_stacklevel,
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
        elif backend == "jax":
            warnings.warn(
                "MRVI model is being loaded with JAX backend",
                DeprecationWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            return JaxMRVI.load(
                dir_path,
                adata=adata,
                accelerator=accelerator,
                device=device,
                prefix=prefix,
                backup_url=backup_url,
                datamodule=datamodule,
            )
        else:
            raise ValueError(f"Unknown backend '{backend}'. Use 'torch' or 'jax'.")
