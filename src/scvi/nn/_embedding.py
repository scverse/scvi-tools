from __future__ import annotations

from typing import TYPE_CHECKING

import torch
from torch import nn

if TYPE_CHECKING:
    from collections.abc import Callable


def _partial_freeze_hook_factory(freeze: int) -> Callable[[torch.Tensor], torch.Tensor]:
    """Factory for a hook that freezes the first ``freeze`` entries in the gradient.

    Parameters
    ----------
    freeze
        Specifies the number of entries to freeze in the gradient.
    """

    def _partial_freeze_hook(grad: torch.Tensor) -> torch.Tensor:
        grad = grad.clone()
        grad[:freeze] = 0.0
        return grad

    return _partial_freeze_hook


class Embedding(nn.Embedding):
    """``EXPERIMENTAL`` Embedding layer with utility methods for extending.

    Notes
    -----
    Lifecycle: experimental in v1.2.
    """

    @classmethod
    def extend(
        cls,
        embedding: Embedding,
        init: int | list[int],
        freeze_prev: bool = True,
    ) -> Embedding:
        """Factory class method for extending an existing :class:`~scvi.nn.Embedding`.

        Initializes new embeddings with random values or copies from the original embedding.

        Parameters
        ----------
        embedding
            Embedding layer to extend. Not modified in-place as a new instance is returned.
        init
            If an ``int``, specifies the number of new random embeddings to initialize with
            :func:`~torch.nn.init.normal_`. If a ``list[int]``, initializes ``len(init)`` new
            embeddings with the values of the given indices in the original embedding.
        freeze_prev
            If ``True``, gradients for the original embeddings are set to zero.

        Returns
        -------
        New :class:`~scvi.nn.Embedding` instance with :attr:`~scvi.nn.Embedding.num_embeddings`
        increased by the specified ``init`` and with gradients frozen for the first
        ``embedding.num_embeddings`` entries if ``freeze_prev`` is ``True``.
        """
        if isinstance(init, int):
            if init <= 0:
                raise ValueError(f"`init` must be greater than 0, got {init}")
            weight = torch.empty((init, embedding.embedding_dim), device=embedding.weight.device)
            nn.init.normal_(weight)
        elif isinstance(init, list):
            weight = embedding.weight[init].clone()
        else:
            raise TypeError(f"`init` must be an `int` or a `list[int]`, got {type(init)}")

        weight = torch.cat([embedding.weight.clone(), weight], dim=0)
        new_embedding = cls(
            num_embeddings=weight.shape[0],
            embedding_dim=weight.shape[1],
            _weight=weight,
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

    def _load_from_state_dict(self, state_dict: dict[str, torch.Tensor], *args, **kwargs) -> None:
        """Load from a state dict. Overrides initialization parameters with the state dict.

        This is necessary because model constructors will pass in the original parameters, which
        will not match the state dict if the embedding was extended. This method overrides the
        correct attributes based on the state dict.
        """
        key = [key for key in state_dict.keys() if "weight" in key][0]
        self.weight = nn.Parameter(state_dict[key])
        self.num_embeddings, self.embedding_dim = self.weight.shape[0], self.weight.shape[1]
        return super()._load_from_state_dict(state_dict, *args, **kwargs)
