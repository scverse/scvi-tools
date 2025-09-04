from __future__ import annotations

from typing import TYPE_CHECKING

from jax import jit

if TYPE_CHECKING:
    from jax import Array
    from jax.typing import ArrayLike


@jit
def rowwise_max_excluding_diagonal(matrix: ArrayLike) -> Array:
    """Get the rowwise maximum of a matrix excluding the diagonal."""
    import jax.numpy as jnp

    assert matrix.ndim == 2
    num_cols = matrix.shape[1]
    mask = (1 - jnp.eye(num_cols)).astype(bool)
    return (jnp.where(mask, matrix, -jnp.inf)).max(axis=1)


def simple_reciprocal(w: ArrayLike, eps: float = 1e-6) -> Array:
    """Convert distances to similarities via a reciprocal."""
    return 1.0 / (w + eps)
