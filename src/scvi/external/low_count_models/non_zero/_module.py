"""Non-zero SCVI VAE module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS
from scvi.module._vae import VAE
from scvi.module.base import LossOutput
from scvi.utils import unsupported_if_adata_minified

if TYPE_CHECKING:
    from torch.distributions import Distribution


class NonZeroVAE(VAE):
    """VAE that only computes reconstruction loss on non-zero observations.

    This variant of the standard VAE addresses the observation that zeros dominate
    training gradients in standard SCVI. By masking zeros from the reconstruction
    loss, the model focuses on accurately predicting non-zero expression values.

    Parameters
    ----------
    n_input
        Number of input features.
    normalize_by_nonzero
        If ``True``, normalize reconstruction loss by the number of non-zero
        observations per cell. This maintains a per-gene scale comparable to
        standard SCVI. If ``False``, simply sum the masked log probabilities.
    **kwargs
        Additional keyword arguments passed to :class:`~scvi.module.VAE`.

    Notes
    -----
    When ``normalize_by_nonzero=True``, both the reconstruction loss AND the KL
    divergence are normalized by the number of non-zero genes per cell. This
    maintains a consistent KL-to-reconstruction ratio across cells with different
    sparsity levels, preventing posterior collapse that would otherwise occur
    when KL (summed over latent dims) dominates the per-gene reconstruction loss.
    """

    def __init__(
        self,
        n_input: int,
        *args,
        normalize_by_nonzero: bool = True,
        **kwargs,
    ):
        super().__init__(n_input, *args, **kwargs)
        self.normalize_by_nonzero = normalize_by_nonzero

    @unsupported_if_adata_minified
    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, Distribution | None],
        kl_weight: torch.Tensor | float = 1.0,
    ) -> LossOutput:
        """Compute the loss with reconstruction only on non-zero observations."""
        from torch.distributions import kl_divergence

        x = tensors[REGISTRY_KEYS.X_KEY]

        kl_divergence_z = kl_divergence(
            inference_outputs[MODULE_KEYS.QZ_KEY],
            generative_outputs[MODULE_KEYS.PZ_KEY],
        ).sum(dim=-1)

        if not self.use_observed_lib_size:
            kl_divergence_l = kl_divergence(
                inference_outputs[MODULE_KEYS.QL_KEY],
                generative_outputs[MODULE_KEYS.PL_KEY],
            ).sum(dim=1)
        else:
            kl_divergence_l = torch.zeros_like(kl_divergence_z)

        px = generative_outputs[MODULE_KEYS.PX_KEY]
        log_prob = px.log_prob(x)  # (batch_size, n_genes)

        nonzero_mask = x > 0
        masked_log_prob = log_prob * nonzero_mask.float()

        if self.normalize_by_nonzero:
            n_nonzero = nonzero_mask.sum(dim=-1).clamp(min=1).float()
            reconst_loss = -masked_log_prob.sum(dim=-1) / n_nonzero
            kl_divergence_z_scaled = kl_divergence_z / n_nonzero
            kl_divergence_l_scaled = kl_divergence_l / n_nonzero
        else:
            reconst_loss = -masked_log_prob.sum(dim=-1)
            kl_divergence_z_scaled = kl_divergence_z
            kl_divergence_l_scaled = kl_divergence_l

        weighted_kl_local = kl_weight * kl_divergence_z_scaled + kl_divergence_l_scaled
        loss = torch.mean(reconst_loss + weighted_kl_local)

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local={
                MODULE_KEYS.KL_L_KEY: kl_divergence_l,
                MODULE_KEYS.KL_Z_KEY: kl_divergence_z,
            },
        )
