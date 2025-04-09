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
