import torch
from pyro.distributions import BetaBinomial as BetaBinomialDistribution
from torch.distributions import constraints
from torch.distributions.utils import broadcast_all

from ._constraints import open_interval, optional_constraint


class BetaBinomial(BetaBinomialDistribution):
    r"""``BETA`` Beta binomial distribution.

    One of the following parameterizations must be provided:

    1. (``alpha``, ``beta``, ``total_counts``) where alpha and beta are the shape parameters of
    the beta distribution and total_counts is the number of trials.

    2. (``mu``, ``gamma``, ``total_counts``), which is the one used by methylVI. These
    parameters respectively control the mean and dispersion of the distribution.

    In the (``mu``, ``gamma``) parameterization, samples from the beta-binomial are generated
    as follows:

    1. :math:`p_i \sim \textrm{Beta}(\mu, \gamma)`
    2. :math:`y_i \sim \textrm{Ber}(p_i)`
    3. :math:`y = \sum_{i}y_i`

    Parameters
    ----------
    total_count
        Number of trials. Must be a non-negative integer.
    alpha
        As in :class:`~pyro.distributions.BetaBinomial`,
        serves as the first shape parameterization of the beta
        distribution. Must be greater than ``0``.
    beta
        As in :class:`~pyro.distributions.BetaBinomial`,
        serves as the second shape parameterization of the beta
        distribution. Must be greater than ``0``.
    mu
        Mean of the distribution. Must be within the interval ``(0, 1)``.
    gamma
        Dispersion. Must be within the interval ``(0, 1)``.
    validate_args
        Raise ``ValueError`` if arguments do not match the constraints.
    eps
        Numerical stability constant. (See Notes)

    Notes
    -----
    Under the hood we use :class:`~pyro.distributions.BetaBinomial` to implement
    the Beta-Binomial distribution. Thus, when the user specifies a (``mu``, ``gamma``)
    parameterization, we must convert to the (``alpha``, ``beta``) parameterization
    used by the underlying Pyro distribution class. During this process, numerical
    stability issues sometimes cause ``alpha`` or ``beta`` to be equal to (exactly) zero.
    This is not allowed (``alpha`` and ``beta`` must be strictly greater than 0), so we clamp these
    values to be greater than a small constant ``eps``.
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
        alpha: torch.Tensor | None = None,
        beta: torch.Tensor | None = None,
        mu: torch.Tensor | None = None,
        gamma: torch.Tensor | None = None,
        validate_args: bool = False,
        eps: float = 1e-8,
    ):
        self._eps = eps

        using_param_1 = alpha is not None and beta is not None
        using_param_2 = mu is not None and gamma is not None

        using_both_params = using_param_1 and using_param_2
        using_neither_param = (not using_param_1) and (not using_param_2)

        if using_both_params or using_neither_param:
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

    def __repr__(self) -> str:
        param_names = [k for k, _ in self.arg_constraints.items() if k in self.__dict__]
        args_string = ", ".join(
            [
                f"{p}: "
                f"{self.__dict__[p] if self.__dict__[p].numel() == 1 else self.__dict__[p].size()}"
                for p in param_names
                if self.__dict__[p] is not None
            ]
        )
        return self.__class__.__name__ + "(" + args_string + ")"
