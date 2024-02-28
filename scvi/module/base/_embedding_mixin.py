import torch
from torch.nn import ModuleDict

from scvi.module.base._decorators import auto_move_data
from scvi.nn import Embedding


class EmbeddingModuleMixin:
    @property
    def embeddings_dict(self) -> ModuleDict:
        """Dictionary of embeddings."""
        if not hasattr(self, "_embeddings_dict"):
            self._embeddings_dict = ModuleDict()
        return self._embeddings_dict

    def add_embedding(self, key: str, embedding: Embedding, overwrite: bool = False) -> None:
        """Add an embedding to the module."""
        if key in self.embeddings_dict and not overwrite:
            raise KeyError(f"Embedding {key} already exists.")
        self.embeddings_dict[key] = embedding

    def remove_embedding(self, key: str) -> None:
        """Remove an embedding from the module."""
        if key not in self.embeddings_dict:
            raise KeyError(f"Embedding {key} not found.")
        del self.embeddings_dict[key]

    def get_embedding(self, key: str) -> Embedding:
        """Get an embedding from the module."""
        if key not in self.embeddings_dict:
            raise KeyError(f"Embedding {key} not found.")
        return self.embeddings_dict[key]

    def init_embedding(
        self,
        key: str,
        num_embeddings: int,
        embedding_dim: int = 5,
        **kwargs,
    ) -> None:
        """Initialize an embedding in the module."""
        self.add_embedding(key, Embedding(num_embeddings, embedding_dim, **kwargs))

    @auto_move_data
    def compute_embedding(self, indices: torch.Tensor, key: str) -> torch.Tensor:
        """Forward pass for an embedding."""
        return self.get_embedding(key)(indices)
