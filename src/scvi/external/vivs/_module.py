from __future__ import annotations

from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.distributions import Bernoulli, Normal

from scvi import REGISTRY_KEYS
from scvi.external.vivs._constants import VIVS_REGISTRY_KEYS
from scvi.module import VAE
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

if TYPE_CHECKING:
    from typing import Literal


class ImportanceScoreNet(nn.Module):
    """Predicts ``Y`` from log-CPM-normalized ``X``.

    Ported from ``ImportanceScorer``/``ImportanceScorerLinear`` in VIVS's original JAX
    implementation (Boyeau et al. 2024, :cite:p:`BoyeauVIVS24`). The per-cell, per-response
    negative log-likelihood (``all_loss``) is the conditional-randomization-test statistic.
    """

    def __init__(
        self,
        n_input: int,
        n_responses: int,
        n_hidden: int = 128,
        dropout_rate: float = 0.0,
        loss_type: Literal["mse", "binary"] = "mse",
        linear: bool = False,
    ):
        super().__init__()
        self.loss_type = loss_type
        self.linear = linear
        # momentum=0.01 here matches flax's momentum=0.99 (flax's `momentum` is the
        # weight kept on the *old* running stat; torch's is the weight given to the *new* batch).
        if linear:
            self.norm1 = nn.BatchNorm1d(n_input, momentum=0.01, eps=1e-3)
            self.dropout1 = nn.Dropout(0.0)
            self.dense_out = nn.Linear(n_input, n_responses)
        else:
            self.dense1 = nn.Linear(n_input, n_hidden)
            self.norm1 = nn.BatchNorm1d(n_hidden, momentum=0.01, eps=1e-3)
            self.dropout1 = nn.Dropout(dropout_rate)
            self.dense_out = nn.Linear(n_hidden, n_responses)
        # Fixed (non-learned) log-std for the "mse" (Normal) likelihood, matching the
        # original's `self.log_std = 0.0` plain attribute (not a trained flax param).
        self.register_buffer("log_std", torch.zeros(()))

    def forward(self, x: torch.Tensor, y: torch.Tensor) -> dict[str, torch.Tensor]:
        if self.linear:
            h = self.norm1(x)
            h = self.dropout1(h)
            h = self.dense_out(h)
        else:
            h = self.dense1(x)
            h = self.norm1(h)
            h = torch.nn.functional.leaky_relu(h)
            h = self.dropout1(h)
            h = self.dense_out(h)

        if self.loss_type == "mse":
            all_loss = -Normal(h, torch.exp(self.log_std)).log_prob(y)
        else:
            all_loss = -Bernoulli(logits=h).log_prob(y)
        return {"h": h, "loss": all_loss.mean(), "all_loss": all_loss}


class VIVSModule(BaseModuleClass):
    """Compound module for VIVS: a VAE over ``X`` plus an importance-score net for ``Y|X``.

    ``self.x_module`` (a :class:`~scvi.module.VAE`) is used only to sample conditional
    "knockoff" replacements of ``X`` for the CRT; ``self.xy_module`` (an
    ``ImportanceScoreNet``) predicts ``Y`` from ``X`` and its
    negative log-likelihood is the CRT test statistic. The two are trained sequentially,
    not jointly (see :class:`~scvi.external.VIVS`), so ``self._phase`` switches which
    component ``loss()`` optimizes.
    """

    def __init__(
        self,
        n_input: int,
        n_responses: int,
        n_batch: int = 0,
        x_model_kwargs: dict | None = None,
        xy_model_kwargs: dict | None = None,
        xy_linear: bool = False,
        xy_include_batch_in_input: bool = False,
        x_module: BaseModuleClass | None = None,
    ):
        super().__init__()
        x_model_kwargs = x_model_kwargs or {}
        xy_model_kwargs = xy_model_kwargs or {}

        self.x_module_is_pretrained = x_module is not None
        if x_module is not None:
            self.x_module = x_module
            self.x_module.requires_grad_(False)
            self.x_module.eval()
        else:
            self.x_module = VAE(n_input=n_input, n_batch=n_batch, **x_model_kwargs)

        self.n_batch = n_batch
        self.xy_include_batch_in_input = xy_include_batch_in_input
        xy_n_input = n_input + n_batch if xy_include_batch_in_input else n_input
        self.xy_module = ImportanceScoreNet(
            n_input=xy_n_input,
            n_responses=n_responses,
            linear=xy_linear,
            **xy_model_kwargs,
        )
        self._phase: str = "x"

    def _get_inference_input(self, tensors):
        return self.x_module._get_inference_input(tensors)

    def _get_generative_input(self, tensors, inference_outputs):
        return self.x_module._get_generative_input(tensors, inference_outputs)

    @auto_move_data
    def inference(self, *args, **kwargs):
        return self.x_module.inference(*args, **kwargs)

    @auto_move_data
    def generative(self, *args, **kwargs):
        return self.x_module.generative(*args, **kwargs)

    @staticmethod
    def normalize_log_cpm(x: torch.Tensor) -> torch.Tensor:
        """Log-CPM normalization shared by every CRT statistic computation.

        Clamps the per-cell library size away from 0: group-level knockoff substitution
        (in ``get_hier_importance``) can zero out a sparse cell's entire count vector when
        the knockoff draw for a whole gene-cluster happens to sample all zeros, which would
        otherwise produce a 0/0 NaN here and crash the downstream ``Normal(...).log_prob``.
        """
        library = x.sum(dim=-1, keepdim=True).clamp_min(1e-6)
        return torch.log1p(1e6 * x / library)

    def xy_input(self, x: torch.Tensor, batch_index: torch.Tensor) -> torch.Tensor:
        """Build the importance-score net's input: normalized X, optionally + one-hot batch."""
        x_norm = self.normalize_log_cpm(x)
        if self.xy_include_batch_in_input:
            batch_oh = torch.nn.functional.one_hot(
                batch_index.squeeze(-1).long(), self.n_batch
            ).float()
            return torch.cat([x_norm, batch_oh], dim=-1)
        return x_norm

    def loss(self, tensors, inference_outputs, generative_outputs, kl_weight: float = 1.0):
        if self._phase == "x":
            return self.x_module.loss(
                tensors, inference_outputs, generative_outputs, kl_weight=kl_weight
            )

        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[VIVS_REGISTRY_KEYS.Y_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        xy_out = self.xy_module(self.xy_input(x, batch_index), y)
        zeros = torch.zeros_like(xy_out["loss"])
        # xy_out["loss"] is already reduced to a 0-dim scalar (all_loss.mean()), so
        # LossOutput can't infer n_obs_minibatch from reconstruction_loss.shape[0];
        # pass it explicitly instead.
        return LossOutput(
            loss=xy_out["loss"],
            reconstruction_loss=xy_out["loss"],
            kl_local=zeros,
            n_obs_minibatch=x.shape[0],
        )
