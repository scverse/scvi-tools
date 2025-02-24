from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.module.base import EmbeddingModuleMixin

if TYPE_CHECKING:
    import numpy as np
    from anndata import AnnData


class EmbeddingMixin:
    """Mixin class for computing covariate embeddings of a model.

    Must be used with a module that inherits from :class:`~scvi.module.base.EmbeddingModuleMixin`.

    """

    @torch.inference_mode()
    def get_batch_representation(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        batch_size: int | None = None,
        key: str = REGISTRY_KEYS.BATCH_KEY,
        return_mean: bool = True,
    ) -> np.ndarray:
        """Get the batch representation for a given set of indices.

        Parameters
        ----------
        adata
            AnnData object to use.
        indices
            Indices to get the batch representation for.
        batch_size
            Minibatch size for computing the batch representation.
        key
            Setup key to compute the batch representation for.
        return_mean
            Return the mean of the batch representation. Or sample from it.
        """
        if not isinstance(self.module, EmbeddingModuleMixin):
            raise ValueError("The current `module` must inherit from `EmbeddingModuleMixin`.")
        if key not in self.module.embeddings_dim:
            raise ValueError(f"Embedding {key} not found. Enable it during model setup.")

        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        tensors = [
            self.module.compute_embedding(key, tensors[key], return_mean=return_mean)
            for tensors in dataloader
        ]
        return torch.cat(tensors).detach().cpu().numpy()
