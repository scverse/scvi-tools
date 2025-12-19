from __future__ import annotations

import warnings

import torch
import torch.nn.functional as F
from torch.distributions import Distribution, Normal, constraints
from torch.distributions.utils import broadcast_all, lazy_property, logits_to_probs, probs_to_logits

from scvi import settings

from ._constraints import optional_constraint


class LogNormal(Distribution):
    r"""Log-Normal distribution.

    Standard log-normal distribution where log(X) follows a Normal distribution.
    Support is (0, ∞) - this distribution cannot model exact zeros.
    Use Log1pNormal or ZeroInflatedLogNormal for data with zeros.

    In the (`mu`, `sigma`) parameterization, samples from the log-normal are generated as
    follows:

    1. :math:`z \sim \textrm{Normal}(\mu, \sigma)`
    2. :math:`x = \exp(z)`

    The probability density function is:

    .. math::

        f(x; \mu, \sigma) = \frac{1}{x \sigma \sqrt{2\pi}}
        \exp\left(-\frac{(\ln x - \mu)^2}{2\sigma^2}\right)

    Parameters
    ----------
    mu
        Mean of the normal distribution in log space.
    sigma
        Standard deviation of the normal distribution in log space.
    scale
        Normalized mean expression of the distribution.
        This optional parameter is not used in any computations but allows storing
        normalization expression levels.
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    arg_constraints = {
        "mu": optional_constraint(constraints.real),
        "sigma": optional_constraint(constraints.greater_than(0)),
        "scale": optional_constraint(constraints.greater_than_eq(0)),
    }
    support = constraints.positive

    def __init__(
        self,
        mu: torch.Tensor,
        sigma: torch.Tensor,
        scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ):
        self.mu, self.sigma = broadcast_all(mu, sigma)
        self.scale = scale
        super().__init__(validate_args=validate_args)

    @property
    def mean(self) -> torch.Tensor:
        """Mean of the distribution."""
        # E[X] = exp(μ + σ²/2)
        return torch.exp(self.mu + 0.5 * self.sigma**2)

    @property
    def variance(self) -> torch.Tensor:
        """Variance of the distribution."""
        # Var[X] = (exp(σ²) - 1) * exp(2μ + σ²)
        sigma_sq = self.sigma**2
        return (torch.exp(sigma_sq) - 1) * torch.exp(2 * self.mu + sigma_sq)

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: torch.Size | tuple | None = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
        # Sample from Normal in log space, then transform back
        normal_dist = Normal(self.mu, self.sigma)
        log_samples = normal_dist.sample(sample_shape)
        # Transform: exp(log_samples)
        return torch.exp(log_samples)

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """Log probability."""
        if self._validate_args:
            try:
                self._validate_sample(value)
            except ValueError:
                warnings.warn(
                    "The value argument must be within the support of the distribution",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )

        # Transform to log space: y = log(x)
        y = torch.log(value)

        # Log probability in transformed space
        log_prob_normal = (
            -0.5 * torch.log(2 * torch.tensor(torch.pi, device=value.device))
            - torch.log(self.sigma)
            - 0.5 * ((y - self.mu) / self.sigma) ** 2
        )

        # Jacobian correction: d/dx log(x) = 1/x
        log_jacobian = -torch.log(value)

        return log_prob_normal + log_jacobian

    def get_normalized(self, key) -> torch.Tensor:
        """Get normalized values."""
        if key == "mu":
            return self.mu
        elif key == "scale":
            return self.scale
        else:
            raise ValueError(f"normalized key {key} not recognized")

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


class Log1pNormal(Distribution):
    r"""Log1p-Normal distribution.

    Distribution where log(X + 1) follows a Normal distribution.
    This allows modeling data with zeros, as the support is [0, ∞).

    In the (`mu`, `sigma`) parameterization, samples from the log1p-normal are generated as
    follows:

    1. :math:`z \sim \textrm{Normal}(\mu, \sigma)`
    2. :math:`x = \exp(z) - 1`

    The probability density function is:

    .. math::

        f(x; \mu, \sigma) = \frac{1}{(x+1) \sigma \sqrt{2\pi}}
        \exp\left(-\frac{(\ln(x+1) - \mu)^2}{2\sigma^2}\right)

    Parameters
    ----------
    mu
        Mean of the normal distribution in log1p space.
    sigma
        Standard deviation of the normal distribution in log1p space.
    scale
        Normalized mean expression of the distribution.
        This optional parameter is not used in any computations but allows storing
        normalization expression levels.
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    arg_constraints = {
        "mu": optional_constraint(constraints.real),
        "sigma": optional_constraint(constraints.greater_than(0)),
        "scale": optional_constraint(constraints.greater_than_eq(0)),
    }
    support = constraints.nonnegative

    def __init__(
        self,
        mu: torch.Tensor,
        sigma: torch.Tensor,
        scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ):
        self.mu, self.sigma = broadcast_all(mu, sigma)
        self.scale = scale
        super().__init__(validate_args=validate_args)

    @property
    def mean(self) -> torch.Tensor:
        """Mean of the distribution."""
        # E[X] = exp(μ + σ²/2) - 1
        return torch.exp(self.mu + 0.5 * self.sigma**2) - 1

    @property
    def variance(self) -> torch.Tensor:
        """Variance of the distribution."""
        # Var[X] = (exp(σ²) - 1) * exp(2μ + σ²)
        sigma_sq = self.sigma**2
        return (torch.exp(sigma_sq) - 1) * torch.exp(2 * self.mu + sigma_sq)

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: torch.Size | tuple | None = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
        # Sample from Normal in log space, then transform back
        normal_dist = Normal(self.mu, self.sigma)
        log_samples = normal_dist.sample(sample_shape)
        # Transform: expm1(x) = exp(x) - 1, inverse of log1p
        # Clamp to ensure non-negative values (expm1 can be negative when log_samples < 0)
        # TODO: verify if clamp is necessary
        return torch.clamp(torch.expm1(log_samples), min=0.0)

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """Log probability with log1p transformation."""
        if self._validate_args:
            try:
                self._validate_sample(value)
            except ValueError:
                warnings.warn(
                    "The value argument must be within the support of the distribution",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )

        # Transform to log space: y = log(x + 1)
        y = torch.log1p(value)

        # Log probability in transformed space
        log_prob_normal = (
            -0.5 * torch.log(2 * torch.tensor(torch.pi, device=value.device))
            - torch.log(self.sigma)
            - 0.5 * ((y - self.mu) / self.sigma) ** 2
        )

        # Jacobian correction: d/dx log(x + 1) = 1/(x + 1)
        log_jacobian = -torch.log1p(value)

        return log_prob_normal + log_jacobian

    def get_normalized(self, key) -> torch.Tensor:
        """Get normalized values."""
        if key == "mu":
            return self.mu
        elif key == "scale":
            return self.scale
        else:
            raise ValueError(f"normalized key {key} not recognized")

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


class ZeroInflatedLogNormal(LogNormal):
    r"""Zero-inflated log-normal distribution.

    A mixture distribution of a point mass at zero and a log-normal distribution.
    This is a mixed distribution with discrete support at zero and continuous support over (0, ∞).

    In the (`mu`, `sigma`, `zi_logits`) parameterization, samples are generated as follows:

    1. :math:`\pi = \textrm{sigmoid}(\texttt{zi\_logits})`
    2. :math:`b \sim \textrm{Bernoulli}(\pi)`
    3. If :math:`b = 1`: :math:`x = 0`
    4. If :math:`b = 0`: :math:`z \sim \textrm{Normal}(\mu, \sigma)`, :math:`x = \exp(z)`

    The probability mass/density function (mixed discrete-continuous) is:

    .. math::

        p(x; \mu, \sigma, \pi) = \begin{cases}
        \pi & \text{if } x = 0 \\
        (1 - \pi) \cdot f_{LN}(x) & \text{if } x > 0
        \end{cases}

    where :math:`f_{LN}(x; \mu, \sigma) = \frac{1}{x \sigma \sqrt{2\pi}}
    \exp\left(-\frac{(\ln x - \mu)^2}{2\sigma^2}\right)` is the log-normal density.

    Note: The point mass at :math:`x = 0` has probability :math:`\pi`, and the continuous
    log-normal component (for :math:`x > 0`) has probability :math:`(1-\pi)`.

    Parameters
    ----------
    mu
        Mean of the normal distribution in log space.
    sigma
        Standard deviation of the normal distribution in log space.
    zi_logits
        Logits scale of zero inflation probability.
    scale
        Normalized mean expression of the distribution.
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    arg_constraints = {
        "mu": optional_constraint(constraints.real),
        "sigma": optional_constraint(constraints.greater_than(0)),
        "zi_logits": optional_constraint(constraints.real),
        "scale": optional_constraint(constraints.greater_than_eq(0)),
    }
    support = constraints.nonnegative

    def __init__(
        self,
        mu: torch.Tensor,
        sigma: torch.Tensor,
        zi_logits: torch.Tensor,
        scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ):
        super().__init__(
            mu=mu,
            sigma=sigma,
            scale=scale,
            validate_args=validate_args,
        )
        self.zi_logits, self.mu, self.sigma = broadcast_all(zi_logits, self.mu, self.sigma)

    @property
    def mean(self) -> torch.Tensor:
        """Mean of the zero-inflated distribution."""
        pi = self.zi_probs
        ln_mean = super().mean
        return (1 - pi) * ln_mean

    @property
    def variance(self) -> torch.Tensor:
        """Variance of the zero-inflated distribution."""
        pi = self.zi_probs
        ln_mean = super().mean
        ln_var = super().variance
        return (1 - pi) * (ln_var + pi * ln_mean**2)

    @lazy_property
    def zi_logits(self) -> torch.Tensor:
        """ZI logits."""
        return probs_to_logits(self.zi_probs, is_binary=True)

    @lazy_property
    def zi_probs(self) -> torch.Tensor:
        """ZI probabilities."""
        return logits_to_probs(self.zi_logits, is_binary=True)

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: torch.Size | tuple | None = None,
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
                stacklevel=settings.warnings_stacklevel,
            )

        eps = 1e-8

        # Zero-inflation mixture:
        # For x = 0: log P(X=0) = log(pi)
        # For x > 0: log P(X=x) = log(1-pi) + log(p_LN(x))

        # Uses log(sigmoid(x)) = -softplus(-x)
        softplus_zi = F.softplus(-self.zi_logits)

        # Case: x = 0 (point mass)
        # log(pi) = log(sigmoid(zi_logits)) = -softplus(-zi_logits)
        case_zero = -softplus_zi
        mul_case_zero = torch.mul((value < eps).type(torch.float32), case_zero)

        # Case: x > 0 (log-normal component)
        # log(1-pi) = log(sigmoid(-zi_logits)) = -softplus(zi_logits) = -softplus_zi - zi_logits
        # We compute log_prob_ln for all values, but it's only used where value > eps
        # Use safe_value to avoid log(0) in intermediate computation
        safe_value = torch.where(value >= eps, value, torch.ones_like(value))
        y = torch.log(safe_value)

        # Log probability of log-normal component
        log_prob_normal = (
            -0.5 * torch.log(2 * torch.tensor(torch.pi, device=value.device))
            - torch.log(self.sigma)
            - 0.5 * ((y - self.mu) / self.sigma) ** 2
        )

        # Jacobian correction: d/dx log(x) = 1/x
        log_jacobian = -torch.log(safe_value)
        log_prob_ln = log_prob_normal + log_jacobian

        # log((1-pi) * p_LN(x)) = log(1-pi) + log(p_LN(x))
        case_non_zero = -softplus_zi - self.zi_logits + log_prob_ln
        mul_case_non_zero = torch.mul((value > eps).type(torch.float32), case_non_zero)

        res = mul_case_zero + mul_case_non_zero

        return res
