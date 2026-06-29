"""Joint Embedding VAE module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.external.joint_embedding_scvi._utils import (
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
    """VAE with joint embedding loss using binomial thinning and CCO :cite:p:`Svensson26`.

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
        Weight for the CCO loss. Default is 100.0.
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
    """

    def __init__(
        self,
        n_input: int,
        *args,
        joint_embedding_weight: float = 100.0,
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

        # Binomial thinning requires integer count data. The Gaussian likelihood
        # operates on real-valued (e.g. log-normalized) input, which thinning would
        # silently corrupt, so the joint embedding objective is incompatible with it.
        if self.use_joint_embedding and self.gene_likelihood == "normal":
            raise ValueError(
                "JointEmbeddingVAE uses binomial thinning, which requires count data, so "
                "gene_likelihood='normal' is incompatible with use_joint_embedding=True. "
                "Use a count likelihood ('nb', 'zinb', 'poisson') or set "
                "use_joint_embedding=False."
            )

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
        # Get the reconstruction weight
        recon_weight = (
            reconstruction_weight
            if reconstruction_weight is not None
            else self.reconstruction_weight
        )

        # Get base loss components from parent (we'll reweight them)
        base_loss = super().loss(tensors, inference_outputs, generative_outputs, kl_weight)

        # Recompute total loss with reconstruction weight
        # base_loss.loss = mean(reconst_loss + kl_weight * kl_z + kl_l)
        # We want: mean(recon_weight * reconst_loss + kl_weight * kl_z + kl_l)
        # Note: reconstruction_loss is stored as a dict in LossOutput, use dict_sum to get tensor
        reconst_loss = base_loss.dict_sum(base_loss.reconstruction_loss)
        kl_local = base_loss.kl_local
        kl_z = kl_local[MODULE_KEYS.KL_Z_KEY]
        kl_l = kl_local[MODULE_KEYS.KL_L_KEY]

        weighted_recon = recon_weight * reconst_loss
        weighted_kl = kl_weight * kl_z + kl_l
        elbo_loss = torch.mean(weighted_recon + weighted_kl)

        # Joint embedding needs at least two cells to estimate per-dimension
        # statistics (correlation / variance); the default splitter can emit a
        # singleton final minibatch, so fall back to the reweighted ELBO there.
        z_full = inference_outputs[MODULE_KEYS.Z_KEY]
        too_small = z_full.shape[0] < 2

        # If joint embedding is disabled, not training, or the batch is too small,
        # return reweighted ELBO
        if not self.use_joint_embedding or not self.training or too_small:
            return LossOutput(
                loss=elbo_loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_local,
                extra_metrics=base_loss.extra_metrics,
            )

        # Get input data
        x = tensors[REGISTRY_KEYS.X_KEY]

        # Sample per-cell thinning probabilities for realistic library size variation
        p = sample_thinning_probs(x, min_library_size=self.min_library_size)

        # Create thinned view via binomial thinning with per-cell probabilities
        x_thin, _ = binomial_split(x, p=p)

        # Get covariates for inference
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY)
        cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY)

        # Encode thinned view
        inference_thin = self._regular_inference(x_thin, batch_index, cont_covs, cat_covs)
        z_thin = inference_thin[MODULE_KEYS.Z_KEY]

        # Compute CCO loss between thinned and full embeddings (z_full computed above)
        cco_components = cross_correlation_loss(
            z_thin, z_full, self.lambda_off_diag, return_components=True
        )
        cco_loss = cco_components["total"]

        # Compute variance loss for both views (VICReg-style anti-collapse)
        var_loss_full = variance_loss(z_full)
        var_loss_thin = variance_loss(z_thin)
        var_loss = (var_loss_full + var_loss_thin) / 2

        # Apply weights
        je_weight = (
            joint_embedding_weight
            if joint_embedding_weight is not None
            else self.joint_embedding_weight
        )
        var_weight = self.variance_weight

        # Total loss: weighted ELBO + CCO + variance
        total_loss = elbo_loss + je_weight * cco_loss + var_weight * var_loss

        # Combine extra metrics - track CCO components separately
        extra_metrics = {
            **base_loss.extra_metrics,
            "cco_loss": cco_loss.detach(),
            "cco_invariance": cco_components["invariance"].detach(),  # self-supervision
            "cco_redundancy": cco_components["redundancy"].detach(),  # anti-collapse
            "variance_loss": var_loss.detach(),  # VICReg-style anti-collapse
        }

        return LossOutput(
            loss=total_loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local,
            extra_metrics=extra_metrics,
        )
