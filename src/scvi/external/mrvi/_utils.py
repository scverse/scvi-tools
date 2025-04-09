from __future__ import annotations

from typing import TYPE_CHECKING

from jax import jit

from scvi.external.mrvi._types import _ComputeLocalStatisticsRequirements

if TYPE_CHECKING:
    from jax import Array
    from jax.typing import ArrayLike

    from scvi.external.mrvi._types import MRVIReduction


def _parse_local_statistics_requirements(
    reductions: list[MRVIReduction],
    needs_mean_rep: bool = False,
    needs_sampled_rep: bool = False,
    needs_mean_dists: bool = False,
    needs_sampled_dists: bool = False,
    needs_normalized_dists: bool = False,
) -> _ComputeLocalStatisticsRequirements:
    """Get requirements for computing local statistics for a set of reductions."""
    ungrouped_reductions = []
    grouped_reductions = []

    for r in reductions:
        if r.input == "mean_representations":
            needs_mean_rep = True
        elif r.input == "sampled_representations":
            needs_sampled_rep = True
        elif r.input == "mean_distances":
            needs_mean_rep = True
            needs_mean_dists = True
        elif r.input == "sampled_distances":
            needs_sampled_rep = True
            needs_sampled_dists = True
        elif r.input == "normalized_distances":
            needs_sampled_rep = True
            needs_normalized_dists = True
        else:
            raise ValueError(f"Unknown reduction input: {r.input}")

        if r.group_by is None:
            ungrouped_reductions.append(r)
        else:
            grouped_reductions.append(r)

    return _ComputeLocalStatisticsRequirements(
        needs_mean_representations=needs_mean_rep,
        needs_sampled_representations=needs_sampled_rep,
        needs_mean_distances=needs_mean_dists,
        needs_sampled_distances=needs_sampled_dists,
        needs_normalized_distances=needs_normalized_dists,
        grouped_reductions=grouped_reductions,
        ungrouped_reductions=ungrouped_reductions,
    )


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
