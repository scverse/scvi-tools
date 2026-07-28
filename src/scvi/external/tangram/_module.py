from typing import NamedTuple

import torch
from torch import nn

from scvi.module.base import BaseModuleClass, LossOutput


class _TANGRAM_REGISTRY_KEYS_NT(NamedTuple):
    SC_KEY: str = "X"
    SP_KEY: str = "Y"
    DENSITY_KEY: str = "DENSITY"


TANGRAM_REGISTRY_KEYS = _TANGRAM_REGISTRY_KEYS_NT()

EPS = 1e-8


def _cosine_similarity(x: torch.Tensor, y: torch.Tensor, dim: int) -> torch.Tensor:
    numerator = (x * y).sum(dim=dim)
    denominator = torch.clamp(
        torch.linalg.norm(x, dim=dim) * torch.linalg.norm(y, dim=dim),
        min=EPS,
    )
    return numerator / denominator


def _density_criterion(log_y_pred: torch.Tensor, y_true: torch.Tensor) -> torch.Tensor:
    # Kl divergence between the predicted and true distributions
    log_y_true = torch.log(y_true + EPS)
    return (y_true * (log_y_true - log_y_pred)).sum()


class TangramMapper(BaseModuleClass):
    """Torch Tangram Mapper Model."""

    def __init__(
        self,
        n_obs_sc: int,
        n_obs_sp: int,
        lambda_g1: float = 1.0,
        lambda_d: float = 0.0,
        lambda_g2: float = 0.0,
        lambda_r: float = 0.0,
        lambda_count: float = 1.0,
        lambda_f_reg: float = 1.0,
        constrained: bool = False,
        target_count: int | None = None,
    ):
        super().__init__()
        self.n_obs_sc = n_obs_sc
        self.n_obs_sp = n_obs_sp
        self.lambda_g1 = lambda_g1
        self.lambda_d = lambda_d
        self.lambda_g2 = lambda_g2
        self.lambda_r = lambda_r
        self.lambda_count = lambda_count
        self.lambda_f_reg = lambda_f_reg
        self.constrained = constrained
        self.target_count = target_count

        self.mapper_unconstrained = nn.Parameter(torch.randn((self.n_obs_sc, self.n_obs_sp)))

        if self.constrained:
            self.filter_unconstrained = nn.Parameter(torch.randn((self.n_obs_sc, 1)))

    def _get_inference_input(self, tensors: dict[str, torch.Tensor]):
        """Get input for inference."""
        return {}

    def inference(self) -> dict:
        """Run inference model."""
        return {}

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
    ):
        return {}

    def generative(self) -> dict:
        """No generative model here."""
        return {}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
    ):
        """Compute loss."""
        sp = tensors[TANGRAM_REGISTRY_KEYS.SP_KEY]
        sc = tensors[TANGRAM_REGISTRY_KEYS.SC_KEY]
        mapper = torch.softmax(self.mapper_unconstrained, dim=1)

        if self.constrained:
            filter = torch.sigmoid(self.filter_unconstrained)
            mapper_filtered = mapper * filter

        if self.lambda_d > 0:
            density = tensors[TANGRAM_REGISTRY_KEYS.DENSITY_KEY].ravel()
            if self.constrained:
                d_pred = torch.log(mapper_filtered.sum(dim=0) / (filter.sum()))
            else:
                d_pred = torch.log(mapper.sum(dim=0) / mapper.shape[0])
            density_term = self.lambda_d * _density_criterion(d_pred, density)
        else:
            density_term = 0

        if self.constrained:
            sc = sc * filter

        g_pred = mapper.T @ sc

        # Expression term
        if self.lambda_g1 > 0:
            gv_term = self.lambda_g1 * _cosine_similarity(sp, g_pred, dim=0).mean()
        else:
            gv_term = 0
        if self.lambda_g2 > 0:
            vg_term = self.lambda_g1 * _cosine_similarity(sp, g_pred, dim=1).mean()
            vg_term = self.lambda_g2 * vg_term
        else:
            vg_term = 0

        expression_term = gv_term + vg_term

        # Regularization terms
        if self.lambda_r > 0:
            regularizer_term = self.lambda_r * (torch.log(mapper) * mapper).sum()
        else:
            regularizer_term = 0

        if self.lambda_count > 0 and self.constrained:
            if self.target_count is None:
                raise ValueError("target_count must be set if in constrained mode.")
            count_term = self.lambda_count * torch.abs(filter.sum() - self.target_count)
        else:
            count_term = 0

        if self.lambda_f_reg > 0 and self.constrained:
            f_reg_t = filter - torch.square(filter)
            f_reg = self.lambda_f_reg * f_reg_t.sum()
        else:
            f_reg = 0

        # Total loss
        total_loss = -expression_term - regularizer_term + count_term + f_reg
        total_loss = total_loss + density_term

        return LossOutput(
            loss=total_loss,
            n_obs_minibatch=sp.shape[0],
            extra_metrics={
                "expression_term": expression_term,
                "regularizer_term": regularizer_term,
            },
        )

    def sample(self, *args, **kwargs):
        """Not implemented for Tangram."""
        raise NotImplementedError
