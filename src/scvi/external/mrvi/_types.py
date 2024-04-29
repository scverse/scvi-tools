from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal

import xarray as xr

PRNGKey = Any


@dataclass(frozen=True)
class MrVIReduction:
    """
    A dataclass object that represents a single reduction for ``MrVI.compute_local_statistics``.

    Parameters
    ----------
    name
        Name of the reduction. Used as the key for the corresponding DataArray.
    input
        Type of input data.
    fn
        Function that computes the reduction.
    group_by
        Covariate name by which to average the computed statistics by. If ``None``,
        the outputs are left at the per-cell granularity.
    """

    name: str
    input: Literal[
        "mean_representations",
        "mean_distances",
        "sampled_representations",
        "sampled_distances",
        "normalized_distances",
    ]
    fn: callable[[xr.DataArray], xr.DataArray] = lambda x: xr.DataArray(x)
    group_by: str | None = None


@dataclass(frozen=True)
class _ComputeLocalStatisticsRequirements:
    """Utility class for the summarized requirements for ``MrVI.compute_local_statistics``."""

    needs_mean_representations: bool
    needs_mean_distances: bool
    needs_sampled_representations: bool
    needs_sampled_distances: bool
    needs_normalized_distances: bool
    ungrouped_reductions: Iterable[MrVIReduction]
    grouped_reductions: Iterable[MrVIReduction]
