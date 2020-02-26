from typing import Optional, Union, Tuple

import torch
from torch.distributions import NegativeBinomial, constraints
import torch.nn.functional as F
from torch.distributions.utils import (
    broadcast_all,
    probs_to_logits,
    lazy_property,
    logits_to_probs,
)


def from_nb_params2_to_1(mu, theta, eps=1e-6):
    r"""
        :param mu: mean of the NB distribution.
        :param theta: overdispersion
        :param eps: constant used for numerical log stability
        Returns the number of failures until the experiment is stopped
            and the success probability.
    """
    assert (mu is None) == (
        theta is None
    ), "If using the mu/theta NB parameterization, both parameters must be specified"
    logits = (mu + eps).log() - (theta + eps).log()
    total_count = theta
    return total_count, logits


def from_nb_params1_to_2(total_count, logits):
    """
        :param total_count: Number of failures until the experiment is stopped.
        :param logits: success logits
        Returns the mean and overdispersion of the NB distribution.
    """
    theta = total_count
    mu = logits.exp() * theta
    return mu, theta


class NB(NegativeBinomial):
    r"""
        Negative Binomial distribution using two parameterizations:
        - The first one corresponds to the parameterization NB(`total_count`, `probs`)
            where `total_count` is the number of failures until the experiment is stopped
            and `probs` the success probability.
        - The (`mu`, `theta`) parameterization is the one used by scVI. These parameters respectively
        control the mean and overdispersion of the distribution.

        `from_nb_params2_to_1` and `from_nb_params1_to_2` provide ways to convert one parameterization to another.

    """

    def __init__(
        self,
        total_count: Optional[Union[torch.Tensor, int, float]] = None,
        probs: Optional[Union[torch.Tensor, float]] = None,
        logits: Optional[Union[torch.Tensor, float]] = None,
        mu: Optional[Union[torch.Tensor, float]] = None,
        theta: Optional[Union[torch.Tensor, float]] = None,
        validate_args=None,
    ):
        self._eps = 1e-8
        if (mu is None) == (total_count is None):
            raise ValueError(
                "Please use one of the two possible parameterizations. Refer to the documentation for more information."
            )
        if mu is not None:
            total_count, logits = from_nb_params2_to_1(mu, theta, eps=self._eps)

        super().__init__(
            total_count, probs=probs, logits=logits, validate_args=validate_args
        )


class ZINB(NegativeBinomial):
    r"""
        Zero Inflated Negative Binomial distribution.
        zi_logits correspond to the zero-inflation logits

        The negative binomial component parameters can follow two two parameterizations:
        - The first one corresponds to the parameterization NB(`total_count`, `probs`)
            where `total_count` is the number of failures until the experiment is stopped
            and `probs` the success probability.
        - The (`mu`, `theta`) parameterization is the one used by scVI. These parameters respectively
        control the mean and overdispersion of the distribution.

        `from_nb_params2_to_1` and `from_nb_params1_to_2` provide ways to convert one parameterization to another.
    """
    arg_constraints = {
        "total_count": constraints.greater_than_eq(0),
        "probs": constraints.half_open_interval(0.0, 1.0),
        "logits": constraints.real,
        "zi_probs": constraints.half_open_interval(0.0, 1.0),
        "zi_logits": constraints.real,
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        total_count: Optional[Union[torch.Tensor, int, float]] = None,
        probs: Optional[Union[torch.Tensor, float]] = None,
        logits: Optional[Union[torch.Tensor, float]] = None,
        mu: Optional[Union[torch.Tensor, float]] = None,
        theta: Optional[Union[torch.Tensor, float]] = None,
        zi_logits: Optional[Union[torch.Tensor, float]] = None,
        validate_args=None,
    ):
        self._eps = 1e-8  # Used for numerical stability

        # Step 1: Harmonizing parameterizations
        if (mu is None) == (total_count is None):
            raise ValueError(
                "Please use one of the two possible parameterizations. Refer to the documentation for more information."
            )
        if mu is not None:
            total_count, logits = from_nb_params2_to_1(mu, theta, eps=self._eps)

        # Step 2: Taking care of the distribution parameters
        if (probs is None) == (logits is None):
            raise ValueError(
                "Either `probs` or `logits` must be specified, but not both."
            )
        if probs is not None:
            self.total_count, self.probs, self.zi_logits = broadcast_all(
                total_count, probs, zi_logits
            )
            self.total_count = self.total_count.type_as(self.probs)
        else:
            self.total_count, self.logits, self.zi_logits = broadcast_all(
                total_count, logits, zi_logits
            )
            self.total_count = self.total_count.type_as(self.logits)

        self._param = self.probs if probs is not None else self.logits
        batch_shape = self._param.size()
        super(NegativeBinomial, self).__init__(batch_shape, validate_args=validate_args)

    @lazy_property
    def zi_logits(self) -> torch.Tensor:
        return probs_to_logits(self.zi_probs, is_binary=True)

    @lazy_property
    def zi_probs(self) -> torch.Tensor:
        return logits_to_probs(self.zi_logits, is_binary=True)

    def sample(
        self, sample_shape: Union[torch.Size, Tuple] = torch.Size()
    ) -> torch.Tensor:
        with torch.no_grad():
            rate = self._gamma.sample(sample_shape=sample_shape)
            samp = torch.poisson(rate)
            is_zero = torch.rand_like(samp) <= self.zi_probs
            samp[is_zero] = 0.0
            return samp

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        log_nb = super().log_prob(value)

        # motivation: log(sigmoid(x)) = -softplus(-x)
        log_p_zi = -F.softplus(-self.zi_logits)
        log_p_not_zi = -F.softplus(self.zi_logits)

        # case 1: obs is zero
        log_obs_nozi = log_nb + log_p_not_zi
        contrib_zero = (
            torch.logsumexp(torch.stack([log_p_zi, log_obs_nozi], dim=0), dim=0)
            * (value <= self._eps).float()
        )

        # case 2: non zero observation ==> necessarily comes from NB
        contrib_nonzero = log_obs_nozi * (value >= self._eps).float()
        return contrib_zero + contrib_nonzero
