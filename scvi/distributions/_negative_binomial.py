import warnings
from typing import Optional, Tuple, Union

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
import torch
import torch.nn.functional as F
from numpyro.distributions import constraints as numpyro_constraints
from numpyro.distributions.util import promote_shapes, validate_sample
from torch.distributions import Distribution, Gamma
from torch.distributions import Poisson as PoissonTorch
from torch.distributions import constraints
from torch.distributions.utils import (
    broadcast_all,
    lazy_property,
    logits_to_probs,
    probs_to_logits,
)


def log_zinb_positive(
    x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor, pi: torch.Tensor, eps=1e-8
):
    """
    Log likelihood (scalar) of a minibatch according to a zinb model.

    Parameters
    ----------
    x
        Data
    mu
        mean of the negative binomial (has to be positive support) (shape: minibatch x vars)
    theta
        inverse dispersion parameter (has to be positive support) (shape: minibatch x vars)
    pi
        logit of the dropout parameter (real support) (shape: minibatch x vars)
    eps
        numerical stability constant

    Notes
    -----
    We parametrize the bernoulli using the logits, hence the softplus functions appearing.
    """
    # theta is the dispersion rate. If .ndimension() == 1, it is shared for all cells (regardless of batch or labels)
    if theta.ndimension() == 1:
        theta = theta.view(
            1, theta.size(0)
        )  # In this case, we reshape theta for broadcasting

    # Uses log(sigmoid(x)) = -softplus(-x)
    softplus_pi = F.softplus(-pi)
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


def log_nb_positive(
    x: Union[torch.Tensor, jnp.ndarray],
    mu: Union[torch.Tensor, jnp.ndarray],
    theta: Union[torch.Tensor, jnp.ndarray],
    eps: float = 1e-8,
    log_fn: callable = torch.log,
    lgamma_fn: callable = torch.lgamma,
):
    """
    Log likelihood (scalar) of a minibatch according to a nb model.

    Parameters
    ----------
    x
        data
    mu
        mean of the negative binomial (has to be positive support) (shape: minibatch x vars)
    theta
        inverse dispersion parameter (has to be positive support) (shape: minibatch x vars)
    eps
        numerical stability constant
    """
    log = log_fn
    lgamma = lgamma_fn
    log_theta_mu_eps = log(theta + mu + eps)
    res = (
        theta * (log(theta + eps) - log_theta_mu_eps)
        + x * (log(mu + eps) - log_theta_mu_eps)
        + lgamma(x + theta)
        - lgamma(theta)
        - lgamma(x + 1)
    )

    return res


def log_mixture_nb(
    x: torch.Tensor,
    mu_1: torch.Tensor,
    mu_2: torch.Tensor,
    theta_1: torch.Tensor,
    theta_2: torch.Tensor,
    pi_logits: torch.Tensor,
    eps=1e-8,
):
    """
    Log likelihood (scalar) of a minibatch according to a mixture nb model.

    pi_logits is the probability (logits) to be in the first component.
    For totalVI, the first component should be background.

    Parameters
    ----------
    x
        Observed data
    mu_1
        Mean of the first negative binomial component (has to be positive support) (shape: minibatch x features)
    mu_2
        Mean of the second negative binomial (has to be positive support) (shape: minibatch x features)
    theta_1
        First inverse dispersion parameter (has to be positive support) (shape: minibatch x features)
    theta_2
        Second inverse dispersion parameter (has to be positive support) (shape: minibatch x features)
        If None, assume one shared inverse dispersion parameter.
    pi_logits
        Probability of belonging to mixture component 1 (logits scale)
    eps
        Numerical stability constant
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

    logsumexp = torch.logsumexp(torch.stack((log_nb_1, log_nb_2 - pi_logits)), dim=0)
    softplus_pi = F.softplus(-pi_logits)

    log_mixture_nb = logsumexp - softplus_pi

    return log_mixture_nb


def _convert_mean_disp_to_counts_logits(mu, theta, eps=1e-6):
    r"""
    NB parameterizations conversion.

    Parameters
    ----------
    mu
        mean of the NB distribution.
    theta
        inverse overdispersion.
    eps
        constant used for numerical log stability. (Default value = 1e-6)

    Returns
    -------
    type
        the number of failures until the experiment is stopped
        and the success probability.
    """
    if not (mu is None) == (theta is None):
        raise ValueError(
            "If using the mu/theta NB parameterization, both parameters must be specified"
        )
    logits = (mu + eps).log() - (theta + eps).log()
    total_count = theta
    return total_count, logits


def _convert_counts_logits_to_mean_disp(total_count, logits):
    """
    NB parameterizations conversion.

    Parameters
    ----------
    total_count
        Number of failures until the experiment is stopped.
    logits
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


class Poisson(PoissonTorch):
    """
    Poisson distribution.

    Parameters
    ----------
    rate
        rate of the Poisson distribution.
    validate_args
        whether to validate input.
    scale
        Normalized mean expression of the distribution.
        This optional parameter is not used in any computations, but allows to store
        normalization expression levels.

    """

    def __init__(
        self,
        rate: torch.Tensor,
        validate_args: Optional[bool] = None,
        scale: Optional[torch.Tensor] = None,
    ):
        super().__init__(rate=rate, validate_args=validate_args)
        self.scale = scale


class NegativeBinomial(Distribution):
    r"""
    Negative binomial distribution.

    One of the following parameterizations must be provided:

    (1), (`total_count`, `probs`) where `total_count` is the number of failures until
    the experiment is stopped and `probs` the success probability. (2), (`mu`, `theta`)
    parameterization, which is the one used by scvi-tools. These parameters respectively
    control the mean and inverse dispersion of the distribution.

    In the (`mu`, `theta`) parameterization, samples from the negative binomial are generated as follows:

    1. :math:`w \sim \textrm{Gamma}(\underbrace{\theta}_{\text{shape}}, \underbrace{\theta/\mu}_{\text{rate}})`
    2. :math:`x \sim \textrm{Poisson}(w)`

    Parameters
    ----------
    total_count
        Number of failures until the experiment is stopped.
    probs
        The success probability.
    mu
        Mean of the distribution.
    theta
        Inverse dispersion.
    scale
        Normalized mean expression of the distribution.
    validate_args
        Raise ValueError if arguments do not match constraints
    """

    arg_constraints = {
        "mu": constraints.greater_than_eq(0),
        "theta": constraints.greater_than_eq(0),
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        total_count: Optional[torch.Tensor] = None,
        probs: Optional[torch.Tensor] = None,
        logits: Optional[torch.Tensor] = None,
        mu: Optional[torch.Tensor] = None,
        theta: Optional[torch.Tensor] = None,
        scale: Optional[torch.Tensor] = None,
        validate_args: bool = False,
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
        self.scale = scale
        super().__init__(validate_args=validate_args)

    @property
    def mean(self):  # noqa: D102
        return self.mu

    @property
    def variance(self):  # noqa: D102
        return self.mean + (self.mean**2) / self.theta

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: Optional[Union[torch.Size, Tuple]] = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
        gamma_d = self._gamma()
        p_means = gamma_d.sample(sample_shape)

        # Clamping as distributions objects can have buggy behaviors when
        # their parameters are too high
        l_train = torch.clamp(p_means, max=1e8)
        counts = PoissonTorch(
            l_train
        ).sample()  # Shape : (n_samples, n_cells_batch, n_vars)
        return counts

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:  # noqa: D102
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
    r"""
    Zero-inflated negative binomial distribution.

    One of the following parameterizations must be provided:

    (1), (`total_count`, `probs`) where `total_count` is the number of failures until
    the experiment is stopped and `probs` the success probability. (2), (`mu`, `theta`)
    parameterization, which is the one used by scvi-tools. These parameters respectively
    control the mean and inverse dispersion of the distribution.

    In the (`mu`, `theta`) parameterization, samples from the negative binomial are generated as follows:

    1. :math:`w \sim \textrm{Gamma}(\underbrace{\theta}_{\text{shape}}, \underbrace{\theta/\mu}_{\text{rate}})`
    2. :math:`x \sim \textrm{Poisson}(w)`

    Parameters
    ----------
    total_count
        Number of failures until the experiment is stopped.
    probs
        The success probability.
    mu
        Mean of the distribution.
    theta
        Inverse dispersion.
    zi_logits
        Logits scale of zero inflation probability.
    scale
        Normalized mean expression of the distribution.
    validate_args
        Raise ValueError if arguments do not match constraints
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
        total_count: Optional[torch.Tensor] = None,
        probs: Optional[torch.Tensor] = None,
        logits: Optional[torch.Tensor] = None,
        mu: Optional[torch.Tensor] = None,
        theta: Optional[torch.Tensor] = None,
        zi_logits: Optional[torch.Tensor] = None,
        scale: Optional[torch.Tensor] = None,
        validate_args: bool = False,
    ):

        super().__init__(
            total_count=total_count,
            probs=probs,
            logits=logits,
            mu=mu,
            theta=theta,
            scale=scale,
            validate_args=validate_args,
        )
        self.zi_logits, self.mu, self.theta = broadcast_all(
            zi_logits, self.mu, self.theta
        )

    @property
    def mean(self):  # noqa: D102
        pi = self.zi_probs
        return (1 - pi) * self.mu

    @property
    def variance(self):  # noqa: D102
        raise NotImplementedError

    @lazy_property
    def zi_logits(self) -> torch.Tensor:
        """ZI logits."""
        return probs_to_logits(self.zi_probs, is_binary=True)

    @lazy_property
    def zi_probs(self) -> torch.Tensor:  # noqa: D102
        return logits_to_probs(self.zi_logits, is_binary=True)

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: Optional[Union[torch.Size, Tuple]] = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
        samp = super().sample(sample_shape=sample_shape)
        is_zero = torch.rand_like(samp) <= self.zi_probs
        samp_ = torch.where(is_zero, torch.zeros_like(samp), samp)
        return samp_

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """Log probability."""
        try:
            self._validate_sample(value)
        except ValueError:
            warnings.warn(
                "The value argument must be within the support of the distribution",
                UserWarning,
            )
        return log_zinb_positive(value, self.mu, self.theta, self.zi_logits, eps=1e-08)


class NegativeBinomialMixture(Distribution):
    """
    Negative binomial mixture distribution.

    See :class:`~scvi.distributions.NegativeBinomial` for further description
    of parameters.

    Parameters
    ----------
    mu1
        Mean of the component 1 distribution.
    mu2
        Mean of the component 2 distribution.
    theta1
        Inverse dispersion for component 1.
    mixture_logits
        Logits scale probability of belonging to component 1.
    theta2
        Inverse dispersion for component 1. If `None`, assumed to be equal to `theta1`.
    validate_args
        Raise ValueError if arguments do not match constraints
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
        validate_args: bool = False,
    ):

        (
            self.mu1,
            self.theta1,
            self.mu2,
            self.mixture_logits,
        ) = broadcast_all(mu1, theta1, mu2, mixture_logits)

        super().__init__(validate_args=validate_args)

        if theta2 is not None:
            self.theta2 = broadcast_all(mu1, theta2)
        else:
            self.theta2 = None

    @property
    def mean(self):  # noqa: D102
        pi = self.mixture_probs
        return pi * self.mu1 + (1 - pi) * self.mu2

    @lazy_property
    def mixture_probs(self) -> torch.Tensor:  # noqa: D102
        return logits_to_probs(self.mixture_logits, is_binary=True)

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: Optional[Union[torch.Size, Tuple]] = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
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
        counts = PoissonTorch(
            l_train
        ).sample()  # Shape : (n_samples, n_cells_batch, n_features)
        return counts

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """Log probability."""
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


class JaxNegativeBinomialMeanDisp(dist.NegativeBinomial2):
    """Negative binomial parameterized by mean and inverse dispersion."""

    arg_constraints = {
        "mean": numpyro_constraints.positive,
        "inverse_dispersion": numpyro_constraints.positive,
    }
    support = numpyro_constraints.nonnegative_integer

    def __init__(
        self,
        mean: jnp.ndarray,
        inverse_dispersion: jnp.ndarray,
        validate_args: Optional[bool] = None,
        eps: float = 1e-8,
    ):
        self._inverse_dispersion, self._mean = promote_shapes(inverse_dispersion, mean)
        self._eps = eps
        super().__init__(mean, inverse_dispersion, validate_args=validate_args)

    @property
    def mean(self):  # noqa: D102
        return self._mean

    @property
    def inverse_dispersion(self):  # noqa: D102
        return self._inverse_dispersion

    @validate_sample
    def log_prob(self, value):
        """Log probability."""
        # theta is inverse_dispersion
        theta = self._inverse_dispersion
        mu = self._mean
        eps = self._eps
        return log_nb_positive(
            value,
            mu,
            theta,
            eps=eps,
            log_fn=jnp.log,
            lgamma_fn=jax.scipy.special.gammaln,
        )
