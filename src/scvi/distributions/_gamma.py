from __future__ import annotations

import warnings

import torch
import torch.nn.functional as F
from torch.distributions import Gamma as GammaTorch
from torch.distributions import constraints
from torch.distributions.utils import broadcast_all, lazy_property, logits_to_probs, probs_to_logits

from scvi import settings

from ._constraints import optional_constraint


class Gamma(GammaTorch):
    r"""Gamma distribution.

    Extends PyTorch's Gamma distribution to include optional scale parameter
    for storing normalized expression levels.

    In the (`concentration`, `rate`) parameterization (shape-rate), samples are drawn
    directly from the Gamma distribution:

    .. math::

        x \sim \textrm{Gamma}(\underbrace{\alpha}_{\text{concentration}},
        \underbrace{\beta}_{\text{rate}})

    The probability density function is:

    .. math::

        f(x; \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} e^{-\beta x}

    Parameters
    ----------
    concentration
        Shape parameter (α > 0) of the Gamma distribution.
    rate
        Rate parameter (β > 0) of the Gamma distribution.
        Note: rate = 1/scale in PyTorch's parameterization.
    scale
        Normalized mean expression of the distribution.
        This optional parameter is not used in any computations but allows storing
        normalization expression levels.
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    arg_constraints = {
        "concentration": constraints.greater_than(0),
        "rate": constraints.greater_than(0),
    }
    support = constraints.nonnegative

    def __init__(
        self,
        concentration: torch.Tensor,
        rate: torch.Tensor,
        scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ):
        super().__init__(concentration=concentration, rate=rate, validate_args=validate_args)
        self.scale = scale

    def get_normalized(self, key) -> torch.Tensor:
        """Get normalized values."""
        if key == "concentration":
            return self.concentration
        elif key == "rate":
            return self.rate
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


class ZeroInflatedGamma(Gamma):
    r"""Zero-inflated Gamma distribution.

    A mixture distribution of a point mass at zero and a Gamma distribution.
    This is ideal for continuous positive data with excess zeros.

    In the (`concentration`, `rate`, `zi_logits`) parameterization, samples are generated
    as follows:

    1. :math:`\pi = \textrm{sigmoid}(\texttt{zi\_logits})`
    2. :math:`b \sim \textrm{Bernoulli}(\pi)`
    3. If :math:`b = 1`: :math:`x = 0`
    4. If :math:`b = 0`: :math:`x \sim \textrm{Gamma}(\alpha, \beta)`

    The probability density function is:

    .. math::

        f(x; \alpha, \beta, \pi) = \begin{cases}
        \pi & \text{if } x = 0 \\
        (1 - \pi) \cdot f_{\Gamma}(x) & \text{if } x > 0
        \end{cases}

    where :math:`f_{\Gamma}(x; \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)}
    x^{\alpha-1} e^{-\beta x}` is the Gamma density.

    Parameters
    ----------
    concentration
        Shape parameter (α > 0) of the Gamma distribution.
    rate
        Rate parameter (β > 0) of the Gamma distribution.
    zi_logits
        Logits scale of zero inflation probability.
    scale
        Normalized mean expression of the distribution.
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    # Note: scale is intentionally not in arg_constraints because it's optional
    arg_constraints = {
        "concentration": optional_constraint(constraints.greater_than(0)),
        "rate": optional_constraint(constraints.greater_than(0)),
        "zi_logits": optional_constraint(constraints.real),
    }
    support = constraints.nonnegative

    def __init__(
        self,
        concentration: torch.Tensor,
        rate: torch.Tensor,
        zi_logits: torch.Tensor,
        scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ):
        super().__init__(
            concentration=concentration,
            rate=rate,
            scale=scale,
            validate_args=validate_args,
        )
        self.zi_logits, self.concentration, self.rate = broadcast_all(
            zi_logits, self.concentration, self.rate
        )

    @property
    def mean(self) -> torch.Tensor:
        """Mean of the zero-inflated distribution."""
        pi = self.zi_probs
        gamma_mean = self.concentration / self.rate
        return (1 - pi) * gamma_mean

    @property
    def variance(self) -> torch.Tensor:
        """Variance of the zero-inflated distribution."""
        pi = self.zi_probs
        gamma_mean = self.concentration / self.rate
        gamma_var = self.concentration / (self.rate**2)
        return (1 - pi) * (gamma_var + pi * gamma_mean**2)

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
        # For x > 0: log P(X=x) = log(1-pi) + log(p_Gamma(x))

        # Uses log(sigmoid(x)) = -softplus(-x)
        softplus_zi = F.softplus(-self.zi_logits)

        # Case: x = 0 (point mass)
        # log(pi) = log(sigmoid(zi_logits)) = -softplus(-zi_logits)
        case_zero = -softplus_zi
        mul_case_zero = torch.mul((value < eps).type(torch.float32), case_zero)

        # Case: x > 0 (Gamma component)
        # log(1-pi) = log(sigmoid(-zi_logits)) = -softplus(zi_logits) = -softplus_zi - zi_logits
        # Use safe_value to avoid log(0) in intermediate computation
        safe_value = torch.where(value >= eps, value, torch.ones_like(value))

        # Log probability of Gamma component
        # log p(x; α, β) = (α-1)*log(x) + α*log(β) - β*x - log(Γ(α))
        log_prob_gamma = (
            (self.concentration - 1) * torch.log(safe_value)
            + self.concentration * torch.log(self.rate)
            - self.rate * value
            - torch.lgamma(self.concentration)
        )

        # log((1-pi) * p_Gamma(x)) = log(1-pi) + log(p_Gamma(x))
        case_non_zero = -softplus_zi - self.zi_logits + log_prob_gamma
        mul_case_non_zero = torch.mul((value > eps).type(torch.float32), case_non_zero)

        res = mul_case_zero + mul_case_non_zero

        return res
