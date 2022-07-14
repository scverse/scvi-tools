from typing import Optional

import torch
from torch.distributions import Bernoulli as BernoulliTorch
class Bernoulli(BernoulliTorch):
    """
    Bernoulli distribution.

    Parameters
    ----------
    probs
        probability p for the bernoulli distribution.
    logits
        log-odds of sampling.
    validate_args
        whether to validate input.

    """
    
    def __init__(
        self,
        probs: torch.Tensor = None,
        logits: torch.Tensor = None,
        validate_args: Optional[bool] = None,
    ):
        super().__init__(probs, logits, validate_args)
