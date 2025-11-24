from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from xarray import DataArray

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Literal


@dataclass(frozen=True)
class MRVIReduction:
    """Reduction dataclass for :meth:`~scvi.external.MRVI.compute_local_statistics`.

    Parameters
    ----------
    name
        Name of the reduction. Used as the key for the corresponding :class:`~xarray.DataArray`.
    input
        Type of input data, must be one of the following:

        * ``"mean_representations"``
        * ``"mean_distances"``
        * ``"sampled_representations"``
        * ``"sampled_distances"``
        * ``"normalized_distances"``
    fn
        Function that computes the reduction.
    group_by
        Covariate name by which to average the computed statistics by. If ``None``, the outputs are
        left at the per-cell granularity.
    """

    name: str
    input: Literal[
        "mean_representations",
        "mean_distances",
        "sampled_representations",
        "sampled_distances",
        "normalized_distances",
    ]
    fn: Callable[[DataArray], DataArray] = lambda x: DataArray(x)
    group_by: str | None = None


@dataclass(frozen=True)
class _ComputeLocalStatisticsRequirements:
    """Summarized requirements for :meth:`~scvi.external.MRVI.compute_local_statistics`."""

    needs_mean_representations: bool
    needs_mean_distances: bool
    needs_sampled_representations: bool
    needs_sampled_distances: bool
    needs_normalized_distances: bool
    ungrouped_reductions: Iterable[MRVIReduction]
    grouped_reductions: Iterable[MRVIReduction]


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
