"""Joint Embedding VAE module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.external.low_count_models.utils import (
    binomial_split,
    cross_correlation_loss,
    sample_thinning_probs,
    variance_loss,
)
from scvi.module._constants import MODULE_KEYS
from scvi.module._vae import VAE
from scvi.module.base import LossOutput

if TYPE_CHECKING:
    from torch.distributions import Distribution


class JointEmbeddingVAE(VAE):
    """VAE with joint embedding loss using binomial thinning and CCO.

    This module extends the standard VAE with a cross-correlation objective (CCO)
    loss that encourages the embedding of a thinned view to match the embedding
    of the original data. This promotes robustness to count dropout/noise.

    Thinning probabilities are dynamically sampled per cell to produce target
    library sizes that are log-uniform between min_library_size and the observed
    library size. This matches realistic library size variation in single-cell data.

    Parameters
    ----------
    n_input
        Number of input features.
    joint_embedding_weight
        Weight for the CCO loss. Default is 1.0.
    lambda_off_diag
        Off-diagonal penalty in CCO loss. Default is 0.01.
    min_library_size
        Minimum target library size for thinning. Default is 10.
        Thinned library sizes are sampled log-uniformly between this
        value and the observed library size.
    reconstruction_weight
        Weight for reconstruction loss. Default is 1.0.
        Set to 0.0 for pure self-supervised training with only CCO loss.
    variance_weight
        Weight for variance regularization loss (VICReg-style). Default is 0.0.
        Set to positive value (e.g., 1.0) to prevent dimension collapse
        in self-supervised training.
    use_joint_embedding
        Whether to use joint embedding loss. Default is True.
        Set to False to train as standard SCVI.
    **kwargs
        Additional keyword arguments passed to :class:`~scvi.module.VAE`.

    See Also
    --------
    :class:`~scvi.module.VAE`
    :class:`~scvi.external.JointEmbeddingSCVI`
    """

    def __init__(
        self,
        n_input: int,
        *args,
        joint_embedding_weight: float = 1.0,
        lambda_off_diag: float = 0.01,
        min_library_size: float = 10.0,
        reconstruction_weight: float = 1.0,
        variance_weight: float = 0.0,
        use_joint_embedding: bool = True,
        **kwargs,
    ):
        super().__init__(n_input, *args, **kwargs)
        self.joint_embedding_weight = joint_embedding_weight
        self.lambda_off_diag = lambda_off_diag
        self.min_library_size = min_library_size
        self.reconstruction_weight = reconstruction_weight
        self.variance_weight = variance_weight
        self.use_joint_embedding = use_joint_embedding

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, Distribution | None],
        kl_weight: float = 1.0,
        joint_embedding_weight: float | None = None,
        reconstruction_weight: float | None = None,
    ) -> LossOutput:
        """Compute the loss including optional joint embedding CCO loss.

        Parameters
        ----------
        tensors
            Dictionary of input tensors.
        inference_outputs
            Dictionary of inference outputs.
        generative_outputs
            Dictionary of generative outputs.
        kl_weight
            Weight for KL divergence term.
        joint_embedding_weight
            Optional override for joint embedding weight.
        reconstruction_weight
            Optional override for reconstruction weight.
            Set to 0.0 for pure self-supervised training.

        Returns
        -------
        LossOutput
            Loss output with total loss, reconstruction loss, KL terms,
            and CCO loss in extra_metrics.
        """
        recon_weight = (
            reconstruction_weight
            if reconstruction_weight is not None
            else self.reconstruction_weight
        )

        base_loss = super().loss(tensors, inference_outputs, generative_outputs, kl_weight)

        reconst_loss = base_loss.dict_sum(base_loss.reconstruction_loss)
        kl_local = base_loss.kl_local
        kl_z = kl_local[MODULE_KEYS.KL_Z_KEY]
        kl_l = kl_local[MODULE_KEYS.KL_L_KEY]

        weighted_recon = recon_weight * reconst_loss
        weighted_kl = kl_weight * kl_z + kl_l
        elbo_loss = torch.mean(weighted_recon + weighted_kl)

        if not self.use_joint_embedding or not self.training:
            return LossOutput(
                loss=elbo_loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_local,
                extra_metrics=base_loss.extra_metrics,
            )

        x = tensors[REGISTRY_KEYS.X_KEY]
        p = sample_thinning_probs(x, min_library_size=self.min_library_size)
        x_thin, _ = binomial_split(x, p=p)

        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY)
        cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY)

        inference_thin = self._regular_inference(x_thin, batch_index, cont_covs, cat_covs)
        z_thin = inference_thin[MODULE_KEYS.Z_KEY]
        z_full = inference_outputs[MODULE_KEYS.Z_KEY]

        cco_components = cross_correlation_loss(
            z_thin, z_full, self.lambda_off_diag, return_components=True
        )
        cco_loss = cco_components["total"]

        var_loss_full = variance_loss(z_full)
        var_loss_thin = variance_loss(z_thin)
        var_loss = (var_loss_full + var_loss_thin) / 2

        je_weight = (
            joint_embedding_weight
            if joint_embedding_weight is not None
            else self.joint_embedding_weight
        )
        var_weight = self.variance_weight

        total_loss = elbo_loss + je_weight * cco_loss + var_weight * var_loss

        extra_metrics = {
            **base_loss.extra_metrics,
            "cco_loss": cco_loss.detach(),
            "cco_invariance": cco_components["invariance"].detach(),
            "cco_redundancy": cco_components["redundancy"].detach(),
            "variance_loss": var_loss.detach(),
        }

        return LossOutput(
            loss=total_loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local,
            extra_metrics=extra_metrics,
        )
