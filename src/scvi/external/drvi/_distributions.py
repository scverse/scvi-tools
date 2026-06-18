from __future__ import annotations

import torch
from torch.distributions import Distribution
from torch.nn import functional as F


class LogNegativeBinomial(Distribution):
    r"""Negative binomial distribution parameterized in log space.

    A negative binomial whose mean and dispersion are supplied as ``log_m`` (:math:`\log\mu`) and
    ``log_r`` (:math:`\log\theta`). The log-probability is the standard negative-binomial
    log-probability, but rewritten so that it is evaluated directly from the log-parameters via
    :func:`~torch.nn.functional.softplus` — this avoids ever materializing :math:`\mu` or
    :math:`\theta` and is numerically stable for very small/large means. This is DRVI's "pnb"
    likelihood and is what enables the additive (log-space) split decoder.

    Mathematically equivalent to :class:`scvi.distributions.NegativeBinomial` with
    ``mu = exp(log_m)``, ``theta = exp(log_r)``.

    Parameters
    ----------
    log_m
        Log of the mean :math:`\log\mu`.
    log_r
        Log of the inverse-dispersion :math:`\log\theta` (``theta = exp(log_r)``).
    log_scale
        Optional log of the library-size-independent normalized mean, exposed as ``scale`` for
        :class:`~scvi.model.base.RNASeqMixin` (e.g. ``log_softmax`` of the decoder output).
    eps
        Small constant for numerical stability.
    """

    arg_constraints = {
        "log_m": torch.distributions.constraints.real,
        "log_r": torch.distributions.constraints.real,
    }
    support = torch.distributions.constraints.nonnegative_integer

    def __init__(
        self,
        log_m: torch.Tensor,
        log_r: torch.Tensor,
        log_scale: torch.Tensor | None = None,
        eps: float = 1e-8,
        validate_args: bool = False,
    ) -> None:
        self._eps = eps
        self.log_m = log_m
        self.log_r = log_r
        # mean-space parameters, exposed as attributes for RNASeqMixin compatibility. log_prob uses
        # log_m / log_r directly (full log-space numerical stability); ``mu`` is a property whose
        # setter keeps ``log_m`` in sync, because RNASeqMixin reassigns ``px.mu`` before calling
        # log_prob (importance weighting / differential expression).
        self._mu = torch.exp(log_m)
        self.theta = torch.exp(log_r)
        self.scale = torch.exp(log_scale) if log_scale is not None else self._mu
        super().__init__(validate_args=validate_args)

    @property
    def mu(self) -> torch.Tensor:
        return self._mu

    @mu.setter
    def mu(self, value: torch.Tensor) -> None:
        self._mu = value
        self.log_m = torch.log(value + self._eps)

    @property
    def mean(self) -> torch.Tensor:
        return self._mu

    @property
    def variance(self) -> torch.Tensor:
        return self._mu + self._mu**2 / self.theta

    def get_normalized(self, key: str) -> torch.Tensor:
        """Return a named mean-space parameter (RNASeqMixin contract)."""
        if key in ("mu", "rate"):
            return self.mu
        elif key == "scale":
            return self.scale
        elif key == "theta":
            return self.theta
        raise ValueError(f"normalized key {key} not recognized")

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """Negative-binomial log-probability evaluated from the log-parameters."""
        log_m, log_r, eps = self.log_m, self.log_r, self._eps
        r = torch.exp(log_r)
        # log C(value + r - 1, value)
        choice = (
            torch.lgamma(value + r + eps) - torch.lgamma(value + 1 + eps) - torch.lgamma(r + eps)
        )
        # value * log(p),  p = mu / (mu + theta);  log(p) = -softplus(log_r - log_m)
        log_pow_k = -value * F.softplus(log_r - log_m + eps)
        # r * log(1 - p),  1 - p = theta / (mu + theta);  log(1 - p) = -softplus(log_m - log_r)
        log_pow_r = -r * F.softplus(log_m - log_r + eps)
        return choice + log_pow_k + log_pow_r

    def sample(self, sample_shape: torch.Size = torch.Size()) -> torch.Tensor:
        """Sample via the Gamma-Poisson mixture (same as a negative binomial).

        Parameterized from the log-parameters to stay consistent with :meth:`log_prob`:
        ``concentration = theta = exp(log_r)`` and ``rate = theta / mu = exp(log_r - log_m)``.
        """
        with torch.no_grad():
            concentration = torch.exp(self.log_r).clamp(min=self._eps)
            rate = torch.exp(self.log_r - self.log_m).clamp(min=self._eps)
            gamma = torch.distributions.Gamma(concentration=concentration, rate=rate)
            samples = gamma.sample(sample_shape)
            return torch.distributions.Poisson(samples.clamp(max=1e8)).sample()
