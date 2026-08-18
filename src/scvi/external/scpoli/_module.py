"""scPoli module: VAE with epoch-averaged cell-type prototype embeddings."""

from __future__ import annotations

from typing import TYPE_CHECKING

import torch

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    import numpy as np
    from torch.distributions import Distribution

from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS
from scvi.module._vae import VAE
from scvi.module.base import LossOutput


class ScPoliVAE(VAE):
    r"""VAE with non-parametric prototype embeddings for cell-type anchoring.

    Extends :class:`~scvi.module.VAE` with learnable cell-type prototypes
    in latent space, following :cite:p:`Lotfollahi23`.

    **Labeled prototypes** are maintained as plain tensors (no gradient flow)
    and recomputed at the end of every prototype-training epoch as the
    per-cell-type mean of latent representations:

    .. math::

        p_k \\leftarrow \\frac{1}{|C_k|}\\sum_{i \\in C_k} z_i

    **Unlabeled prototypes** are initialized by clustering the latent space
    (KMeans or Leiden) and subsequently updated at the end of each epoch via a
    dedicated Adam optimizer that minimizes the mean nearest-centroid distance
    (k-means style objective).  These positions are fixed during each mini-batch
    forward pass; only the encoder is trained toward them.

    Prototype computation is gated by ``_prototypes_initialized`` /
    ``_unlabeled_prototypes_initialized`` flags so that the ELBO-only
    pre-training phase incurs zero overhead.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_labels
        Number of *known* cell-type labels (the unlabeled category must be
        excluded; its integer code is ``>= n_labels``).
    eta
        Weight on the labeled prototype loss term.
    unlabeled_weight
        Weight on the unlabeled prototype loss.  ``0`` (default) disables the
        unlabeled loss entirely; a typical value used in scArches is ``0``.
    **vae_kwargs
        All remaining keyword arguments are forwarded to :class:`~scvi.module.VAE`.
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 2,
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        dropout_rate: float = 0.0,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        encode_covariates: bool = True,
        deeply_inject_covariates: bool = False,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_size_factor_key: bool = False,
        use_observed_lib_size: bool = True,
        library_log_means: np.ndarray | None = None,
        library_log_vars: np.ndarray | None = None,
        var_activation: Callable[[torch.Tensor], torch.Tensor] = None,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        batch_embedding_kwargs: dict | None = None,
        eta: float = 0.0,
        unlabeled_weight: float = 0.01,
    ):
        super().__init__(
            n_input=n_input,
            n_batch=n_batch,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            encode_covariates=encode_covariates,
            deeply_inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            use_size_factor_key=use_size_factor_key,
            use_observed_lib_size=use_observed_lib_size,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            var_activation=var_activation,
            extra_encoder_kwargs=extra_encoder_kwargs,
            extra_decoder_kwargs=extra_decoder_kwargs,
            # batch_representation="embedding" activates VAE's init_embedding()
            # path; must be set explicitly alongside batch_embedding_kwargs.
            batch_representation="embedding" if batch_embedding_kwargs is not None else "one-hot",
            batch_embedding_kwargs=batch_embedding_kwargs,
        )
        self.n_prototypes = n_labels  # one prototype per known cell type
        self._n_latent = n_latent  # stored for dynamic buffer allocation
        self.eta = eta
        self.unlabeled_weight = unlabeled_weight

        # Labeled prototype buffers — registered unconditionally so that
        # _prototypes_initialized is always accessible (even when n_labels=0).
        # Shape (0, n_latent) when there are no labeled cell types.
        self.register_buffer(
            "prototypes_labeled",
            torch.zeros(self.n_prototypes, n_latent),
        )
        self.register_buffer(
            "_prototypes_initialized",
            torch.zeros(1, dtype=torch.bool),
        )

        # Unlabeled prototype buffers — empty until clustering is done
        # Shape (0, n_latent) signals "not yet allocated".
        # set_unlabeled_prototypes() replaces this with the real tensor.
        self.register_buffer(
            "prototypes_unlabeled",
            torch.zeros(0, n_latent),
        )
        self.register_buffer(
            "_unlabeled_prototypes_initialized",
            torch.zeros(1, dtype=torch.bool),
        )

    def _load_from_state_dict(
        self,
        state_dict,
        prefix,
        local_metadata,
        strict,
        missing_keys,
        unexpected_keys,
        error_msgs,
    ):
        """Custom state dict loading to handle dynamic-size prototypes_unlabeled buffer.

        The prototypes_unlabeled buffer can have variable shape [n_clusters, n_latent]
        depending on clustering results. When loading a saved model, we need to resize
        the buffer to match the saved shape before loading the state dict.
        """
        # Check if prototypes_unlabeled is in the state dict and has a different size
        unlabeled_key = prefix + "prototypes_unlabeled"
        if unlabeled_key in state_dict:
            saved_shape = state_dict[unlabeled_key].shape
            current_shape = self.prototypes_unlabeled.shape
            # If shapes differ, re-register the buffer with the correct shape
            if saved_shape != current_shape:
                self.register_buffer(
                    "prototypes_unlabeled",
                    torch.zeros(saved_shape),
                )

        # Call parent implementation
        super()._load_from_state_dict(
            state_dict,
            prefix,
            local_metadata,
            strict,
            missing_keys,
            unexpected_keys,
            error_msgs,
        )

    @torch.no_grad()
    def initialize_prototypes(self, latent: torch.Tensor, labels: torch.Tensor) -> None:
        """Seed each labeled prototype from the mean latent of its cell type.

        Called once at the phase boundary (epoch == ``pretrain_epochs``) before
        any prototype loss is applied, so prototypes start at the true cluster
        centres rather than at zero.

        Parameters
        ----------
        latent
            Tensor of shape ``(n_obs, n_latent)`` — encoder output for training cells.
        labels
            Integer tensor of shape ``(n_obs,)`` with label codes. Cells with
            code ``>= n_prototypes`` are unlabeled and are ignored.
        """
        if self.n_prototypes == 0:
            return
        for k in range(self.n_prototypes):
            mask = labels == k
            if mask.any():
                self.prototypes_labeled[k] = latent[mask].mean(0)
        self._prototypes_initialized.fill_(True)

    @torch.no_grad()
    def update_prototypes(self, latent: torch.Tensor, labels: torch.Tensor) -> None:
        """Recompute each labeled prototype as the mean latent of its cell type.

        Called at the end of every epoch after the prototype phase begins.
        Only prototypes for which labeled training cells exist are updated;
        others keep their previous value (important during query surgery when
        some reference cell types may be absent from the query batch).

        Parameters
        ----------
        latent
            Tensor of shape ``(n_obs, n_latent)``.
        labels
            Integer tensor of shape ``(n_obs,)``.
        """
        if self.n_prototypes == 0:
            return
        for k in range(self.n_prototypes):
            mask = labels == k
            if mask.any():
                self.prototypes_labeled[k] = latent[mask].mean(0)

    def set_unlabeled_prototypes(self, centers: torch.Tensor) -> None:
        """Allocate (or replace) the unlabeled prototype buffer.

        Called from :class:`~scvi.external.scpoli._model.ScPoliPrototypeCallback`
        after clustering the training latent space.  Uses ``register_buffer``
        so the tensor travels with ``model.to(device)`` but is excluded from
        the main optimiser.

        Parameters
        ----------
        centers
            Tensor of shape ``(n_clusters, n_latent)`` with cluster centroids.
        """
        self.register_buffer("prototypes_unlabeled", centers.clone())
        self._unlabeled_prototypes_initialized.fill_(True)

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, Distribution | None],
        kl_weight: float = 1.0,
    ) -> LossOutput:
        """ELBO + labeled and unlabeled prototype losses.

        Prototype losses are only applied after the respective prototype sets
        have been initialised via the training callback.
        """
        # Standard VAE ELBO
        loss_output = super().loss(
            tensors=tensors,
            inference_outputs=inference_outputs,
            generative_outputs=generative_outputs,
            kl_weight=kl_weight,
        )

        proto_loss = torch.tensor(0.0, device=loss_output.loss.device)
        unlabeled_proto_loss = torch.tensor(0.0, device=loss_output.loss.device)

        z = inference_outputs[MODULE_KEYS.Z_KEY]  # (n_cells, n_latent)
        # Labeled prototype loss
        prototypes_ready = (
            self.n_prototypes > 0
            and self._prototypes_initialized.item()
            and REGISTRY_KEYS.LABELS_KEY in tensors
        )
        if prototypes_ready:
            labels = tensors[REGISTRY_KEYS.LABELS_KEY].squeeze(-1).long()
            labeled_mask = labels < self.n_prototypes  # exclude unlabeled cells
            if labeled_mask.any():
                z_labeled = z[labeled_mask]  # (n_lab, n_latent)
                labels_labeled = labels[labeled_mask]  # (n_lab,)

                # Pairwise L2 distances: (n_lab, n_prototypes)
                dists = torch.cdist(z_labeled, self.prototypes_labeled, p=2)

                # Mean L2 distance per cell type (matching scArches formula)
                for k in labels_labeled.unique():
                    idx = (labels_labeled == k).nonzero(as_tuple=False)[:, 0]
                    proto_loss = proto_loss + dists[idx, k].sum() / len(idx)

                proto_loss = proto_loss

        # Unlabeled prototype loss (nearest-centroid / soft k-means)
        # Prototypes are fixed during the forward pass (requires_grad=False);
        # the encoder learns to bring cells closer to the nearest centroid.
        # The centroid positions themselves are updated in the callback via a
        # separate optimizer after each epoch (alternating optimisation).
        unlabeled_ready = (
            self._unlabeled_prototypes_initialized.item()
            and self.unlabeled_weight > 0.0
            and self.prototypes_unlabeled.shape[0] > 0
        )
        if unlabeled_ready:
            # All cells (labeled + unlabeled) contribute to the unlabeled loss
            dists_u = torch.cdist(z, self.prototypes_unlabeled, p=2)
            min_dist, y_hat = torch.min(dists_u, dim=1)
            unique_clusters = y_hat.unique()
            if len(unique_clusters) > 0:
                unlabeled_proto_loss = torch.stack(
                    [min_dist[y_hat == c].mean() for c in unique_clusters]
                ).mean()
                proto_loss = proto_loss + self.unlabeled_weight * unlabeled_proto_loss

        total_loss = loss_output.loss + self.eta * proto_loss
        return LossOutput(
            loss=total_loss,
            reconstruction_loss=loss_output.reconstruction_loss,
            kl_local=loss_output.kl_local,
            extra_metrics={
                "proto_loss": proto_loss,
                "unlabeled_proto_loss": unlabeled_proto_loss,
            },
        )
