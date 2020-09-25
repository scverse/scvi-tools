from typing import Union, Tuple, Optional
import warnings

import torch
import torch.nn.functional as F
from torch.distributions import constraints, Distribution, Gamma, Poisson
from torch.distributions.utils import (
    broadcast_all,
    probs_to_logits,
    lazy_property,
    logits_to_probs,
)


def log_zinb_positive(x, mu, theta, pi, eps=1e-8):
    """Note: All inputs are torch Tensors
    log likelihood (scalar) of a minibatch according to a zinb model.
    Notes:
    We parametrize the bernoulli using the logits, hence the softplus functions appearing

    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    pi: logit of the dropout parameter (real support) (shape: minibatch x genes)
    eps: numerical stability constant

    """

    # theta is the dispersion rate. If .ndimension() == 1, it is shared for all cells (regardless of batch or labels)
    if theta.ndimension() == 1:
        theta = theta.view(
            1, theta.size(0)
        )  # In this case, we reshape theta for broadcasting

    softplus_pi = F.softplus(-pi)  #  uses log(sigmoid(x)) = -softplus(-x)
    log_theta_eps = torch.log(theta + eps)
    log_theta_mu_eps = torch.log(theta + mu + eps)
    pi_theta_log = -pi + theta * (log_theta_eps - log_theta_mu_eps)

    case_zero = F.softplus(pi_theta_log) - softplus_pi
    mul_case_zero = torch.mul((x < eps).type(torch.float32), case_zero)

    case_non_zero = (
        -softplus_pi
        + pi_theta_log
        + x * (torch.log(mu + eps) - log_theta_mu_eps)
        + torch.lgamma(x + theta)
        - torch.lgamma(theta)
        - torch.lgamma(x + 1)
    )
    mul_case_non_zero = torch.mul((x > eps).type(torch.float32), case_non_zero)

    res = mul_case_zero + mul_case_non_zero

    return res


def log_nb_positive(x, mu, theta, eps=1e-8):
    """Note: All inputs should be torch Tensors
    log likelihood (scalar) of a minibatch according to a nb model.

    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    eps: numerical stability constant

    """
    if theta.ndimension() == 1:
        theta = theta.view(
            1, theta.size(0)
        )  # In this case, we reshape theta for broadcasting

    log_theta_mu_eps = torch.log(theta + mu + eps)

    res = (
        theta * (torch.log(theta + eps) - log_theta_mu_eps)
        + x * (torch.log(mu + eps) - log_theta_mu_eps)
        + torch.lgamma(x + theta)
        - torch.lgamma(theta)
        - torch.lgamma(x + 1)
    )

    return res


def log_mixture_nb(x, mu_1, mu_2, theta_1, theta_2, pi, eps=1e-8):
    """Note: All inputs should be torch Tensors
    log likelihood (scalar) of a minibatch according to a mixture nb model.
    pi is the probability to be in the first component.

    For totalVI, the first component should be background.

    Parameters
    ----------
    mu1: mean of the first negative binomial component (has to be positive support) (shape: minibatch x genes)
    theta1: first inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    mu2: mean of the second negative binomial (has to be positive support) (shape: minibatch x genes)
    theta2: second inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
        If None, assume one shared inverse dispersion parameter.
    eps: numerical stability constant


    Returns
    -------

    """
    if theta_2 is not None:
        log_nb_1 = log_nb_positive(x, mu_1, theta_1)
        log_nb_2 = log_nb_positive(x, mu_2, theta_2)
    # this is intended to reduce repeated computations
    else:
        theta = theta_1
        if theta.ndimension() == 1:
            theta = theta.view(
                1, theta.size(0)
            )  # In this case, we reshape theta for broadcasting

        log_theta_mu_1_eps = torch.log(theta + mu_1 + eps)
        log_theta_mu_2_eps = torch.log(theta + mu_2 + eps)
        lgamma_x_theta = torch.lgamma(x + theta)
        lgamma_theta = torch.lgamma(theta)
        lgamma_x_plus_1 = torch.lgamma(x + 1)

        log_nb_1 = (
            theta * (torch.log(theta + eps) - log_theta_mu_1_eps)
            + x * (torch.log(mu_1 + eps) - log_theta_mu_1_eps)
            + lgamma_x_theta
            - lgamma_theta
            - lgamma_x_plus_1
        )
        log_nb_2 = (
            theta * (torch.log(theta + eps) - log_theta_mu_2_eps)
            + x * (torch.log(mu_2 + eps) - log_theta_mu_2_eps)
            + lgamma_x_theta
            - lgamma_theta
            - lgamma_x_plus_1
        )

    logsumexp = torch.logsumexp(torch.stack((log_nb_1, log_nb_2 - pi)), dim=0)
    softplus_pi = F.softplus(-pi)

    log_mixture_nb = logsumexp - softplus_pi

    return log_mixture_nb


def _convert_mean_disp_to_counts_logits(mu, theta, eps=1e-6):
    r"""NB parameterizations conversion

    Parameters
    ----------
    mu :
        mean of the NB distribution.
    theta :
        inverse overdispersion.
    eps :
        constant used for numerical log stability. (Default value = 1e-6)

    Returns
    -------
    type
        the number of failures until the experiment is stopped
        and the success probability.
    """
    assert (mu is None) == (
        theta is None
    ), "If using the mu/theta NB parameterization, both parameters must be specified"
    logits = (mu + eps).log() - (theta + eps).log()
    total_count = theta
    return total_count, logits


def _convert_counts_logits_to_mean_disp(total_count, logits):
    """NB parameterizations conversion

    Parameters
    ----------
    total_count :
        Number of failures until the experiment is stopped.
    logits :
        success logits.

    Returns
    -------
    type
        the mean and inverse overdispersion of the NB distribution.

    """
    theta = total_count
    mu = logits.exp() * theta
    return mu, theta


def _gamma(theta, mu):
    concentration = theta
    rate = theta / mu
    # Important remark: Gamma is parametrized by the rate = 1/scale!
    gamma_d = Gamma(concentration=concentration, rate=rate)
    return gamma_d


class NegativeBinomial(Distribution):
    r"""Negative Binomial(NB) distribution using two parameterizations:

    - (`total_count`, `probs`) where `total_count` is the number of failures
        until the experiment is stopped
        and `probs` the success probability.
    - The (`mu`, `theta`) parameterization is the one used by scVI. These parameters respectively
    control the mean and overdispersion of the distribution.

    `_convert_mean_disp_to_counts_logits` and `_convert_counts_logits_to_mean_disp` provide ways to convert
    one parameterization to another.

    Parameters
    ----------

    Returns
    -------
    """
    arg_constraints = {
        "mu": constraints.greater_than_eq(0),
        "theta": constraints.greater_than_eq(0),
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        total_count: torch.Tensor = None,
        probs: torch.Tensor = None,
        logits: torch.Tensor = None,
        mu: torch.Tensor = None,
        theta: torch.Tensor = None,
        validate_args=True,
    ):
        self._eps = 1e-8
        if (mu is None) == (total_count is None):
            raise ValueError(
                "Please use one of the two possible parameterizations. Refer to the documentation for more information."
            )

        using_param_1 = total_count is not None and (
            logits is not None or probs is not None
        )
        if using_param_1:
            logits = logits if logits is not None else probs_to_logits(probs)
            total_count = total_count.type_as(logits)
            total_count, logits = broadcast_all(total_count, logits)
            mu, theta = _convert_counts_logits_to_mean_disp(total_count, logits)
        else:
            mu, theta = broadcast_all(mu, theta)
        self.mu = mu
        self.theta = theta
        super().__init__(validate_args=validate_args)

    def sample(self, sample_shape=torch.Size()):
        with torch.no_grad():
            gamma_d = self._gamma()
            p_means = gamma_d.sample(sample_shape)

            # Clamping as distributions objects can have buggy behaviors when
            # their parameters are too high
            l_train = torch.clamp(p_means, max=1e8)
            counts = Poisson(
                l_train
            ).sample()  # Shape : (n_samples, n_cells_batch, n_genes)
            return counts

    def log_prob(self, value):
        if self._validate_args:
            try:
                self._validate_sample(value)
            except ValueError:
                warnings.warn(
                    "The value argument must be within the support of the distribution",
                    UserWarning,
                )
        return log_nb_positive(value, mu=self.mu, theta=self.theta, eps=self._eps)

    def _gamma(self):
        return _gamma(self.theta, self.mu)


class ZeroInflatedNegativeBinomial(NegativeBinomial):
    r"""Zero Inflated Negative Binomial distribution.

    zi_logits correspond to the zero-inflation logits
        mu + mu ** 2 / theta
    The negative binomial component parameters can follow two two parameterizations:
    - The first one corresponds to the parameterization NB(`total_count`, `probs`)
        where `total_count` is the number of failures until the experiment is stopped
        and `probs` the success probability.
    - The (`mu`, `theta`) parameterization is the one used by scVI. These parameters respectively
    control the mean and overdispersion of the distribution.

    `_convert_mean_disp_to_counts_logits` and `_convert_counts_logits_to_mean_disp`
    provide ways to convert one parameterization to another.

    Parameters
    ----------

    Returns
    -------
    """
    arg_constraints = {
        "mu": constraints.greater_than_eq(0),
        "theta": constraints.greater_than_eq(0),
        "zi_probs": constraints.half_open_interval(0.0, 1.0),
        "zi_logits": constraints.real,
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        total_count: torch.Tensor = None,
        probs: torch.Tensor = None,
        logits: torch.Tensor = None,
        mu: torch.Tensor = None,
        theta: torch.Tensor = None,
        zi_logits: torch.Tensor = None,
        validate_args=True,
    ):

        super().__init__(
            total_count=total_count,
            probs=probs,
            logits=logits,
            mu=mu,
            theta=theta,
            validate_args=validate_args,
        )
        self.zi_logits, self.mu, self.theta = broadcast_all(
            zi_logits, self.mu, self.theta
        )

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
            samp = super().sample(sample_shape=sample_shape)
            is_zero = torch.rand_like(samp) <= self.zi_probs
            samp[is_zero] = 0.0
            return samp

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        try:
            self._validate_sample(value)
        except ValueError:
            warnings.warn(
                "The value argument must be within the support of the distribution",
                UserWarning,
            )
        return log_zinb_positive(value, self.mu, self.theta, self.zi_logits, eps=1e-08)


class NegativeBinomialMixture(Distribution):
    r"""
    Negative binomial mixture distribution.

    Parameters
    ----------

    Returns
    -------
    """
    arg_constraints = {
        "mu1": constraints.greater_than_eq(0),
        "mu2": constraints.greater_than_eq(0),
        "theta1": constraints.greater_than_eq(0),
        "mixture_probs": constraints.half_open_interval(0.0, 1.0),
        "mixture_logits": constraints.real,
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        mu1: torch.Tensor,
        mu2: torch.Tensor,
        theta1: torch.Tensor,
        mixture_logits: torch.Tensor,
        theta2: Optional[torch.Tensor] = None,
        validate_args=True,
    ):

        (
            self.mu1,
            self.theta1,
            self.mu2,
            self.mixture_logits,
        ) = broadcast_all(mu1, theta1, mu2, mixture_logits)

        if theta2 is not None:
            self.theta2 = broadcast_all(mu1, theta2)
        else:
            self.theta2 = None

        super().__init__(validate_args=validate_args)

    @lazy_property
    def mixture_probs(self) -> torch.Tensor:
        return logits_to_probs(self.mixture_logits, is_binary=True)

    def sample(
        self, sample_shape: Union[torch.Size, Tuple] = torch.Size()
    ) -> torch.Tensor:
        with torch.no_grad():
            pi = self.mixture_probs
            mixing_sample = torch.distributions.Bernoulli(pi).sample()
            mu = self.mu1 * mixing_sample + self.mu2 * (1 - mixing_sample)
            if self.theta2 is None:
                theta = self.theta1
            else:
                theta = self.theta1 * mixing_sample + self.theta2 * (1 - mixing_sample)
            gamma_d = _gamma(mu, theta)
            p_means = gamma_d.sample(sample_shape)

            # Clamping as distributions objects can have buggy behaviors when
            # their parameters are too high
            l_train = torch.clamp(p_means, max=1e8)
            counts = Poisson(
                l_train
            ).sample()  # Shape : (n_samples, n_cells_batch, n_features)
            return counts

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        try:
            self._validate_sample(value)
        except ValueError:
            warnings.warn(
                "The value argument must be within the support of the distribution",
                UserWarning,
            )
        return log_mixture_nb(
            value,
            self.mu1,
            self.mu2,
            self.theta1,
            self.theta2,
            self.mixture_logits,
            eps=1e-08,
        )
