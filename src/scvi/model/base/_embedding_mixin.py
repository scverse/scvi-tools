from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.module.base import EmbeddingModuleMixin

if TYPE_CHECKING:
    import numpy as np
    from anndata import AnnData


class EmbeddingMixin:
    """``EXPERIMENTAL`` Mixin class for initializing and using embeddings in a model.

    Must be used with a module that inherits from :class:`~scvi.module.base.EmbeddingModuleMixin`.

    Notes
    -----
    Lifecycle: experimental in v1.2.
    """

    @torch.inference_mode()
    def get_batch_representation(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        batch_size: int | None = None,
    ) -> np.ndarray:
        """Get the batch representation for a given set of indices."""
        if not isinstance(self.module, EmbeddingModuleMixin):
            raise ValueError("The current `module` must inherit from `EmbeddingModuleMixin`.")

        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        key = REGISTRY_KEYS.BATCH_KEY
        tensors = [self.module.compute_embedding(key, tensors[key]) for tensors in dataloader]
        return torch.cat(tensors).detach().cpu().numpy()
