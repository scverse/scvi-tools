import torch
from torch.distributions import Normal
from torch.nn import ModuleDict

from scvi.module.base._decorators import auto_move_data
from scvi.nn import Embedding


class EmbeddingModuleMixin:
    """Mixin class for initializing and using embeddings in a module."""

    @property
    def embeddings_dict(self) -> ModuleDict:
        """Dictionary of embeddings."""
        if not hasattr(self, "_embeddings_dict"):
            self._embeddings_dict = ModuleDict()
        return self._embeddings_dict

    @property
    def embeddings_dim(self) -> dict:
        """Dictionary of embeddings dimensions."""
        if not hasattr(self, "_embeddings_dict"):
            self._embeddings_dim = {}
        return self._embeddings_dim

    @property
    def variational(self) -> dict:
        """Dictionary of whether embedding is variational."""
        if not hasattr(self, "_variational"):
            self._variational = {}
        return self._variational

    def add_embedding(self, key: str, embedding: Embedding, overwrite: bool = False) -> None:
        """Add an embedding to the module."""
        if key in self.embeddings_dict and not overwrite:
            raise KeyError(f"Embedding {key} already exists.")
        torch.nn.init.zeros_(embedding.weight)
        self.embeddings_dict[key] = embedding

    def remove_embedding(self, key: str) -> None:
        """Remove an embedding from the module."""
        if key not in self.embeddings_dict:
            raise KeyError(f"Embedding {key} not found.")
        del self.embeddings_dict[key]

    def get_embedding(
        self,
        key: str,
    ) -> Embedding:
        """Get an embedding from the module."""
        if key not in self.embeddings_dict:
            raise KeyError(f"Embedding {key} not found.")
        return self.embeddings_dict[key]

    def get_embedding_dim(self, key: str, default_value: str | None = None) -> int:
        """Get the dimension of an embedding."""
        if key not in self.embeddings_dim:
            if default_value is not None:
                return default_value
            else:
                raise KeyError(f"Embedding {key} not found.")
        return self.embeddings_dim[key]

    def get_embedding_variational(self, key: str, default_value: str | None = None) -> bool:
        """Get whether an embedding is variational."""
        if key not in self.variational:
            if default_value is not None:
                return default_value
            else:
                raise KeyError(f"Embedding {key} not found.")
        return self.variational[key]

    def init_embedding(
        self,
        key: str,
        num_embeddings: int,
        embedding_dim: int = 5,
        variational: bool = False,
        **kwargs,
    ) -> None:
        """Initialize an embedding in the module."""
        self.embeddings_dim[key] = embedding_dim
        self.variational[key] = variational

        if variational:
            embedding_dim *= 2
        self.add_embedding(key, Embedding(num_embeddings, embedding_dim, **kwargs))

    @auto_move_data
    def compute_embedding(
        self,
        key: str,
        indices: torch.Tensor,
        return_mean: bool = False,
        return_dist: bool = False,
    ) -> torch.Tensor:
        """Forward pass for an embedding."""
        indices = indices.flatten() if indices.ndim > 1 else indices
        embedding = self.get_embedding(key)(indices)
        if self.get_embedding_variational(key):
            embedding_dim = self.get_embedding_dim(key)
            if return_mean:
                return embedding[:, :embedding_dim]
            dist = Normal(embedding[:, :embedding_dim], torch.exp(embedding[:, embedding_dim:]))
            if return_dist:
                return dist
            else:
                return dist.rsample()
        else:
            return embedding
