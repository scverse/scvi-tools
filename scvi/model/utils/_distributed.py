from __future__ import annotations

import logging

from lightning.pytorch.strategies import DDPStrategy, Strategy

logger = logging.getLogger(__name__)


def use_distributed_sampler(strategy: str | Strategy) -> bool:
    """``EXPERIMENTAL`` Return whether to use a distributed sampler.

    Currently only supports DDP.
    """
    if isinstance(strategy, str):
        # ["ddp", "ddp_spawn", "ddp_find_unused_parameters_true"]
        return "ddp" in strategy
    return isinstance(strategy, DDPStrategy)
