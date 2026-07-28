from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import torch
from tqdm import tqdm

import scvi
from scvi import REGISTRY_KEYS
from scvi.data._utils import _validate_adata_dataloader_input
from scvi.module._constants import MODULE_KEYS

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from anndata import AnnData
    from torch import Tensor

logger = logging.getLogger(__name__)


class GenerativeMixin:
    """Generative helpers backing the DRVI interpretability analyses."""

    @torch.inference_mode()
    def iterate_on_ae_output(
        self,
        adata: AnnData | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        deterministic: bool = False,
    ):
        """Yield ``(inference_outputs, generative_outputs)`` per minibatch.

        Runs the full autoencoder (encoder + decoder) over ``adata`` (or a custom ``dataloader``,
        e.g. from an out-of-core datamodule) via ``module.forward(compute_loss=False)`` (no loss
        computed). With ``deterministic=True`` the bottleneck uses the posterior mean.
        """
        _validate_adata_dataloader_input(self, adata, dataloader)
        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        prev_inspect_mode = self.module.inspect_mode
        prev_fully_deterministic = self.module.fully_deterministic
        self.module.inspect_mode = True
        try:
            if deterministic:
                self.module.fully_deterministic = True
            for tensors in tqdm(dataloader, mininterval=5.0):
                yield self.module.forward(tensors=tensors, compute_loss=False)
        finally:
            self.module.fully_deterministic = prev_fully_deterministic
            self.module.inspect_mode = prev_inspect_mode

    @torch.inference_mode()
    def iterate_on_decoded_latent_samples(
        self,
        z: np.ndarray,
        lib: np.ndarray | None = None,
        batch_values: np.ndarray | None = None,
        cat_values: np.ndarray | None = None,
        cont_values: np.ndarray | None = None,
        batch_size: int = scvi.settings.batch_size,
        map_cat_values: bool = False,
    ):
        """Yield decoder (generative) outputs for given latent samples, batched.

        Parameters
        ----------
        z
            Latent samples, shape ``(n_samples, n_latent)``.
        lib
            Library sizes, shape ``(n_samples,)``. Defaults to ``1e4`` per sample.
        batch_values
            Batch indices, shape ``(n_samples,)``. Defaults to ``0``.
        cat_values
            Categorical covariates, shape ``(n_samples, n_cat_covs)``.
        cont_values
            Continuous covariates, shape ``(n_samples, n_cont_covs)``.
        batch_size
            Minibatch size.
        map_cat_values
            Map categorical / batch values to integer codes via the AnnData manager registries.
        """
        was_training = self.module.training
        prev_inspect_mode = self.module.inspect_mode
        self.module.eval()
        self.module.inspect_mode = True

        if cat_values is not None:
            if cat_values.ndim == 1:
                cat_values = cat_values.reshape(-1, 1)
            if map_cat_values:
                mapped_values = np.zeros_like(cat_values)
                state = self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY)
                for i, (_label, map_keys) in enumerate(state["mappings"].items()):
                    cat_mapping = dict(zip(map_keys, range(len(map_keys)), strict=False))
                    mapped_values[:, i] = np.vectorize(cat_mapping.get)(cat_values[:, i])
                cat_values = mapped_values.astype(np.int32)

        if batch_values is not None:
            batch_values = batch_values.flatten()
            if map_cat_values:
                map_keys = self.adata_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY)[
                    "categorical_mapping"
                ]
                batch_mapping = dict(zip(map_keys, range(len(map_keys)), strict=False))
                batch_values = np.vectorize(batch_mapping.get)(batch_values).astype(np.int32)
            batch_values = batch_values.reshape(-1, 1)

        try:
            for i in np.arange(0, z.shape[0], batch_size):
                sl = np.arange(i, min(i + batch_size, z.shape[0]))
                z_tensor = torch.tensor(z[sl], dtype=torch.float32)
                lib_arr = np.full(sl.shape[0], 1e4) if lib is None else lib[sl]
                library = torch.log(torch.tensor(lib_arr, dtype=torch.float32)).reshape(-1, 1)
                cat_tensor = (
                    torch.tensor(cat_values[sl], dtype=torch.long)
                    if cat_values is not None
                    else None
                )
                cont_tensor = (
                    torch.tensor(cont_values[sl], dtype=torch.float32)
                    if cont_values is not None
                    else None
                )
                batch_tensor = (
                    torch.tensor(batch_values[sl], dtype=torch.long)
                    if batch_values is not None
                    else torch.zeros((sl.shape[0], 1), dtype=torch.long)
                )

                # inject the library through the inference_outputs so we can reuse scvi's
                # VAE._get_generative_input unchanged.
                gen_input = self.module._get_generative_input(
                    tensors={
                        REGISTRY_KEYS.BATCH_KEY: batch_tensor,
                        REGISTRY_KEYS.LABELS_KEY: torch.zeros_like(batch_tensor),
                        REGISTRY_KEYS.CONT_COVS_KEY: cont_tensor,
                        REGISTRY_KEYS.CAT_COVS_KEY: cat_tensor,
                    },
                    inference_outputs={
                        MODULE_KEYS.Z_KEY: z_tensor,
                        MODULE_KEYS.LIBRARY_KEY: library,
                    },
                )
                yield self.module.generative(**gen_input)
        finally:
            self.module.inspect_mode = prev_inspect_mode
            self.module.train(was_training)
