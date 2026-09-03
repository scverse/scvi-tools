from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING

import torch

if TYPE_CHECKING:
    from collections.abc import Callable


@cache
def _mps_supports(op: Callable[[torch.Tensor], torch.Tensor]) -> bool:
    """Whether ``op`` has an MPS kernel in the torch build in use.

    ``aten::_standard_gamma`` gained one in torch 2.12 and ``aten::poisson`` in
    2.14, so the answer depends on the installed version. It is probed rather
    than version-compared because source and nightly builds report versions
    that a comparison reads wrongly.

    Callers must check ``torch.backends.mps.is_available()`` first; this is only
    reached from a distribution whose parameters already live on ``mps``.
    """
    try:
        op(torch.ones(1, device="mps"))
    except (NotImplementedError, RuntimeError):
        return False
    return True


def _needs_cpu_detour(on_mps: bool, op: Callable[[torch.Tensor], torch.Tensor]) -> bool:
    """Whether sampling ``op`` has to be routed through the CPU."""
    return on_mps and not _mps_supports(op)


@cache
def _mps_supports_lgamma_on_noncontiguous() -> bool:
    """Whether ``torch.lgamma`` can run on a non-contiguous (broadcast-expanded) MPS tensor.

    Probed rather than version-compared for the same reason as :func:`_mps_supports`.
    """
    try:
        torch.lgamma(torch.ones(1, 4, device="mps").expand(4, 4))
    except (NotImplementedError, RuntimeError):
        return False
    return True


def subset_distribution(
    my_distribution: torch.distributions.Distribution,
    index: torch.Tensor,
    dim: int = 0,
) -> torch.distributions.Distribution:
    """Utility function to subset the parameter of a Pytorch distribution."""
    return my_distribution.__class__(
        **{
            name: torch.index_select(getattr(my_distribution, name), dim=dim, index=index)
            for name in my_distribution.arg_constraints.keys()
        }
    )


class DistributionConcatenator:
    """Utility class to concatenate Pytorch distributions and move them to cpu.

    All distributions must be of the same type.
    """

    def __init__(self):
        self._params = None
        self.distribution_cls = None

    def store_distribution(self, dist: torch.distributions.Distribution):
        """Add a dictionary of distributions to the concatenator.

        Parameters
        ----------
        dist:
            A Pytorch distribution.
        """
        if self._params is None:
            self._params = {name: [] for name in dist.arg_constraints.keys()}
            self.distribution_cls = dist.__class__
        new_params = {name: getattr(dist, name).cpu() for name in dist.arg_constraints.keys()}
        for param_name, param in new_params.items():
            self._params[param_name].append(param)

    def get_concatenated_distributions(self, axis=0):
        """Returns a concatenated `Distribution` object along the specified axis."""
        concat_params = {key: torch.cat(value, dim=axis) for key, value in self._params.items()}
        return self.distribution_cls(**concat_params)
