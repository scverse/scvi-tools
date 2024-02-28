from __future__ import annotations

import torch
from torch import nn


def _partial_freeze_hook_factory(freeze: int) -> callable:
    """Factory for a hook that partially freezes gradients.

    Parameters
    ----------
    freeze
        Freeze the first ``freeze`` entries in the gradient.
    """

    def _partial_freeze_hook(grad: torch.Tensor) -> torch.Tensor:
        grad = grad.clone()
        grad[:freeze] = 0.0
        return grad

    return _partial_freeze_hook


class Embedding(nn.Embedding):
    """Embedding layer with utility methods for extending."""

    @classmethod
    def extend(
        cls,
        embedding: Embedding,
        init: int | list[int],
        freeze_prev: bool = True,
    ) -> Embedding:
        """Factory method for extending an :class:`~scvi.nn.Embedding`.

        Parameters
        ----------
        embedding
            Embedding layer to extend. Not modified in-place as a new instance is returned.
        init
            If an ``int``, the number of new embeddings to initialize with
            :func:`~torch.nn.init.normal_`. If a ``list[int]``, initializes ``len(init)`` new
            embeddings with the values of the given indices in the original embedding.
        freeze_prev
            If ``True``, gradients for the original embeddings are set to zero.
        """
        old_weight = embedding.weight.clone()  # (E, D)

        if isinstance(init, int):
            if init <= 0:
                raise ValueError(f"`init` must be greater than 0, got {init}")
            n_init = init
            # (F, D)
            new_weight = torch.empty((init, old_weight.shape[1]), device=old_weight.device)
            nn.init.normal_(new_weight)
        elif isinstance(init, list):
            n_init = len(init)
            # (F, D)
            new_weight = old_weight[init]

        new_embedding = cls(
            num_embeddings=embedding.num_embeddings + n_init,
            embedding_dim=embedding.embedding_dim,
            _weight=torch.cat([old_weight, new_weight], dim=0),  # (E + F, D)
            padding_idx=embedding.padding_idx,
            max_norm=embedding.max_norm,
            norm_type=embedding.norm_type,
            scale_grad_by_freq=embedding.scale_grad_by_freq,
            sparse=embedding.sparse,
        )

        if freeze_prev:
            new_embedding.weight.register_hook(
                _partial_freeze_hook_factory(embedding.num_embeddings)
            )

        return new_embedding

    def _load_from_state_dict(self, state_dict: dict[str, torch.Tensor], *args, **kwargs):
        """Load from a state dict. Overrides the initialization parameters with the state dict."""
        weight_tensor = state_dict.get("weight")
        self.weight = nn.Parameter(weight_tensor)
        self.num_embeddings = weight_tensor.shape[0]
        self.embedding_dim = weight_tensor.shape[1]

        return super()._load_from_state_dict(state_dict, *args, **kwargs)
