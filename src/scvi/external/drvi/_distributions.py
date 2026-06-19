from __future__ import annotations

import torch
import torch.nn.functional as F
from torch.distributions import Distribution, constraints
from torch.distributions.utils import broadcast_all

from scvi.distributions._negative_binomial import _gamma, torch_lgamma_mps


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
    validate_args
        Raise ValueError if arguments do not match constraints.
    """

    arg_constraints = {
        "log_m": constraints.real,
        "log_r": constraints.real,
    }
    support = constraints.nonnegative_integer

    def __init__(
        self,
        log_m: torch.Tensor,
        log_r: torch.Tensor,
        log_scale: torch.Tensor | None = None,
        validate_args: bool = False,
    ) -> None:
        self._eps = 1e-8
        self.on_mps = log_m.device.type == "mps"  # TODO: until torch solves the MPS issues
        if log_scale is None:
            self.log_m, self.log_r = broadcast_all(log_m, log_r)
            self.log_scale = None
        else:
            self.log_m, self.log_r, self.log_scale = broadcast_all(log_m, log_r, log_scale)
        super().__init__(validate_args=validate_args)

    # Mean-space parameters are derived lazily from the log-parameters
    # log_prob uses log_m and log_r directly for full log-space numerical stability.
    @property
    def mu(self) -> torch.Tensor:
        return torch.exp(self.log_m)

    @mu.setter
    def mu(self, value: torch.Tensor) -> None:
        # RNASeqMixin reassigns ``px.mu`` before log_prob (importance weighting / DE); keep log_m
        # in sync so log_prob reflects the new mean.
        self.log_m = torch.log(value + self._eps)

    @property
    def theta(self) -> torch.Tensor:
        return torch.exp(self.log_r)

    @property
    def scale(self) -> torch.Tensor:
        return torch.exp(self.log_scale) if self.log_scale is not None else self.mu

    @property
    def mean(self) -> torch.Tensor:
        return self.mu

    @property
    def variance(self) -> torch.Tensor:
        return self.mu + self.mu**2 / self.theta

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
        lgamma = torch_lgamma_mps if self.on_mps else torch.lgamma  # TODO: TORCH MPS FIX
        r = torch.exp(log_r)
        # log C(value + r - 1, value)
        choice = lgamma(value + r + eps) - lgamma(value + 1 + eps) - lgamma(r + eps)
        # value * log(p),  p = mu / (mu + theta);  log(p) = -softplus(log_r - log_m)
        log_pow_k = -value * F.softplus(log_r - log_m + eps)
        # r * log(1 - p),  1 - p = theta / (mu + theta);  log(1 - p) = -softplus(log_m - log_r)
        log_pow_r = -r * F.softplus(log_m - log_r + eps)
        return choice + log_pow_k + log_pow_r

    @torch.inference_mode()
    def sample(self, sample_shape: torch.Size | tuple | None = None) -> torch.Tensor:
        """Sample via the Gamma-Poisson mixture (same as a negative binomial)."""
        sample_shape = sample_shape or torch.Size()
        gamma_d = _gamma(self.theta, self.mu, self.on_mps)  # TODO: TORCH MPS FIX - DONE ON CPU
        p_means = gamma_d.sample(sample_shape)
        # Clamp as the distribution objects can behave badly when their parameters are too high.
        l_train = torch.clamp(p_means, max=1e8)
        counts = (
            torch.distributions.Poisson(l_train).sample().to("mps")
            if self.on_mps  # TODO: NEED TORCH MPS FIX for 'aten::poisson'
            else torch.distributions.Poisson(l_train).sample()
        )
        return counts

    def __repr__(self) -> str:
        param_names = [k for k, _ in self.arg_constraints.items() if k in self.__dict__]
        args_string = ", ".join(
            f"{p}: {v if v.numel() == 1 else v.size()}"
            for p in param_names
            if (v := self.__dict__[p]) is not None
        )
        return self.__class__.__name__ + "(" + args_string + ")"
