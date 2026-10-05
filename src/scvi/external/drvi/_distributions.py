from __future__ import annotations

import torch
import torch.nn.functional as F
from torch.distributions import Distribution, constraints
from torch.distributions.utils import broadcast_all

from scvi.distributions import NegativeBinomial, Normal, Poisson, ZeroInflatedNegativeBinomial
from scvi.distributions._negative_binomial import _gamma, _lgamma_fn
from scvi.distributions._utils import _needs_cpu_detour


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
        self._device = log_m.device
        self.on_mps = self._device.type == "mps"
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
        lgamma = _lgamma_fn(self.on_mps)
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
        gamma_d = _gamma(self.theta, self.mu, self.on_mps)
        p_means = gamma_d.sample(sample_shape)
        # Clamp as the distribution objects can behave badly when their parameters are too high.
        l_train = torch.clamp(p_means, max=1e8)
        if _needs_cpu_detour(self.on_mps, torch.poisson):
            l_train = l_train.to("cpu")
        counts = torch.distributions.Poisson(l_train).sample()
        return counts.to(self._device)

    def __repr__(self) -> str:
        param_names = [k for k, _ in self.arg_constraints.items() if k in self.__dict__]
        args_string = ", ".join(
            f"{p}: {v if v.numel() == 1 else v.size()}"
            for p in param_names
            if (v := self.__dict__[p]) is not None
        )
        return self.__class__.__name__ + "(" + args_string + ")"


def build_gene_likelihood(
    gene_likelihood: str,
    px_scale_logit: torch.Tensor,
    px_r_logit: torch.Tensor,
    px_dropout_logit: torch.Tensor | None = None,
    size_factor: torch.Tensor | None = None,
) -> Distribution:
    """Build the per-gene count/normal distribution from a DRVI decoder's log-space parameters.

    Parameters
    ----------
    gene_likelihood
        One of ``"nb"``, ``"pnb"``, ``"zinb"``, ``"poisson"``, ``"normal"``, ``"normal_unit_var"``.
    px_scale_logit
        Aggregated per-gene log-space scale logits, shape ``(*, n_genes)``.
    px_r_logit
        Per-gene dispersion logit: ``log theta`` for the NB family, ``log variance`` for the
        normal likelihoods.
    px_dropout_logit
        ZINB dropout logits (used only when ``gene_likelihood == "zinb"``).
    size_factor
        Log library size, shape ``(*, 1)``; added to the log-softmax scale to form ``log(mu)``.
        Required for the count likelihoods (``"nb"``, ``"pnb"``, ``"zinb"``, ``"poisson"``); unused
        by ``"normal"`` / ``"normal_unit_var"``.
    """
    if gene_likelihood in ("nb", "pnb", "zinb", "poisson"):
        # log_softmax over genes, then add the (log) library size to get log(mu)
        px_scale_log = px_scale_logit - torch.logsumexp(px_scale_logit, dim=-1, keepdim=True)
        px_rate_log = size_factor + px_scale_log  # size_factor == log(library size)
        if gene_likelihood == "pnb":
            return LogNegativeBinomial(log_m=px_rate_log, log_r=px_r_logit, log_scale=px_scale_log)
        px_scale = torch.exp(px_scale_log)
        px_rate = torch.exp(px_rate_log)
        if gene_likelihood == "nb":
            return NegativeBinomial(mu=px_rate, theta=torch.exp(px_r_logit), scale=px_scale)
        if gene_likelihood == "poisson":
            return Poisson(rate=px_rate, scale=px_scale)
        return ZeroInflatedNegativeBinomial(
            mu=px_rate, theta=torch.exp(px_r_logit), zi_logits=px_dropout_logit, scale=px_scale
        )
    if gene_likelihood == "normal":
        # Gaussian with the mean modeled directly (the raw log-space decoder output, no
        # library/softmax) and per-gene variance modeled in log space.
        var = torch.nan_to_num(torch.exp(px_r_logit), posinf=100.0, neginf=0.0) + 1e-8
        return Normal(px_scale_logit, var.sqrt(), normal_mu=px_scale_logit)
    if gene_likelihood == "normal_unit_var":
        return Normal(px_scale_logit, torch.ones_like(px_scale_logit), normal_mu=px_scale_logit)
    raise ValueError(f"Unknown gene_likelihood: {gene_likelihood}")
