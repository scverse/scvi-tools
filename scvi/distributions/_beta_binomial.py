from typing import Optional

import torch
from pyro.distributions import BetaBinomial as BetaBinomialDistribution
from torch.distributions import constraints
from torch.distributions.utils import broadcast_all

from ._constraints import open_interval, optional_constraint


class BetaBinomial(BetaBinomialDistribution):
    r"""``BETA`` Beta binomial distribution.

    One of the following parameterizations must be provided:

    1. (`alpha`, `beta`, `total_counts`) where alpha and beta are the shape parameters of 
    the beta distribution and total_counts is the number of trials.
    (2), (`mu`, `gamma`, `total_counts`), which is the one used by methylVI. These
    parameters respectively control the mean and dispersion of the distribution.

    In the (`mu`, `gamma`) parameterization, samples from the beta-binomial are generated
    as follows:

    1. :math:`p_i \sim \textrm{Beta}(\mu, \gamma)`
    2. :math:`y_i \sim \textrm{Ber}(p_i)`
    3. :math:`y = \sum_{i}y_i`

    Parameters
    ----------
    total_count
        Number of trials.
    alpha
        First shape parameter of the beta distribution.
    beta
        Second shape parameter of the beta distribution.
    mu
        Mean of the distribution.
    gamma
        Dispersion.
    validate_args
        Raise ValueError if arguments do not match constraints
    eps
        Numerical stability constant
    """

    arg_constraints = {
        "total_count": constraints.nonnegative_integer,
        "alpha": optional_constraint(constraints.greater_than(0)),
        "beta": optional_constraint(constraints.greater_than(0)),
        "mu": optional_constraint(open_interval(0, 1)),
        "gamma": optional_constraint(open_interval(0, 1)),
    }

    support = constraints.nonnegative_integer

    def __init__(
        self,
        total_count: torch.Tensor,
        alpha: Optional[torch.Tensor] = None,
        beta: Optional[torch.Tensor] = None,
        mu: Optional[torch.Tensor] = None,
        gamma: Optional[torch.Tensor] = None,
        validate_args: bool = False,
        eps: float = 1e-8,
    ):
        self._eps = eps

        using_param_1 = alpha is not None and beta is not None
        using_param_2 = mu is not None and gamma is not None

        if (not using_param_1) and (not using_param_2):
            raise ValueError(
                "Please use one of the two possible parameterizations."
                " Refer to the documentation for more information."
            )

        if using_param_1:
            alpha, beta = broadcast_all(alpha, beta)
        else:
            mu, gamma = broadcast_all(mu, gamma)

            alpha = mu * (1 - gamma) / gamma
            beta = (mu - 1) * (gamma - 1) / gamma

            # Clamping to force alpha, beta > 0 due to previously observed
            # numerical stability issues.
            alpha = torch.clamp(alpha, min=self._eps)
            beta = torch.clamp(beta, min=self._eps)

        self.mu = mu
        self.gamma = gamma
        self.alpha = alpha
        self.beta = beta

        super().__init__(
            concentration1=alpha,
            concentration0=beta,
            total_count=total_count,
            validate_args=validate_args,
        )
