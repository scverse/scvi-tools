from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import BaseModuleClass, auto_move_data

if TYPE_CHECKING:
    import numpy as np
    from torch.distributions import Distribution


class VAEC(BaseModuleClass):
    """Conditional Variational auto-encoder model.

    This is an implementation of the CondSCVI model

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches. If ``0``, no batch correction is performed.
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    ct_weight
        Multiplicative weight for cell type specific latent space.
    dropout_rate
        Dropout rate for the encoder and decoder neural network.
    encode_covariates
        If ``True``, covariates are concatenated to gene expression prior to passing through
        the encoder(s). Else, only gene expression is used.
    extra_encoder_kwargs
        Keyword arguments passed into :class:`~scvi.nn.Encoder`.
    extra_decoder_kwargs
        Keyword arguments passed into :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        log_variational: bool = True,
        ct_weight: np.ndarray | None = None,
        dropout_rate: float = 0.05,
        encode_covariates: bool = False,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
    ):
        from scvi.nn import Encoder, FCLayers

        super().__init__()
        self.dispersion = "gene"
        self.n_latent = n_latent
        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.dropout_rate = dropout_rate
        self.log_variational = log_variational
        self.gene_likelihood = "nb"
        self.latent_distribution = "normal"
        self.encode_covariates = encode_covariates
        self.n_batch = n_batch
        self.n_labels = n_labels

        if self.encode_covariates and self.n_batch < 1:
            raise ValueError("`n_batch` must be greater than 0 if `encode_covariates` is `True`.")

        self.px_r = torch.nn.Parameter(torch.randn(n_input))
        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_cat_list=[n_labels] + ([n_batch] if n_batch > 0 and encode_covariates else []),
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            return_dist=True,
            **(extra_encoder_kwargs or {}),
        )

        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=[n_labels] + ([n_batch] if n_batch > 0 else []),
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            **(extra_decoder_kwargs or {}),
        )
        self.px_decoder = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_input), torch.nn.Softplus()
        )

        self.register_buffer(
            "ct_weight",
            (
                torch.ones((self.n_labels,), dtype=torch.float32)
                if ct_weight is None
                else torch.tensor(ct_weight, dtype=torch.float32)
            ),
        )

    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor | None]:
        return {
            MODULE_KEYS.X_KEY: tensors[REGISTRY_KEYS.X_KEY],
            MODULE_KEYS.Y_KEY: tensors[REGISTRY_KEYS.LABELS_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors.get(REGISTRY_KEYS.BATCH_KEY, None),
        }

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution],
    ) -> dict[str, torch.Tensor]:
        return {
            MODULE_KEYS.Z_KEY: inference_outputs[MODULE_KEYS.Z_KEY],
            MODULE_KEYS.LIBRARY_KEY: inference_outputs[MODULE_KEYS.LIBRARY_KEY],
            MODULE_KEYS.Y_KEY: tensors[REGISTRY_KEYS.LABELS_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors.get(REGISTRY_KEYS.BATCH_KEY, None),
        }

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, torch.Tensor | Distribution]:
        """High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        library = x.sum(1).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log1p(x_)

        encoder_input = [x_, y]
        if batch_index is not None and self.encode_covariates:
            encoder_input.append(batch_index)

        qz, z = self.z_encoder(*encoder_input)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand((n_samples, library.size(0), library.size(1)))

        return {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.LIBRARY_KEY: library,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        y: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
    ) -> dict[str, Distribution]:
        """Runs the generative model."""
        from scvi.distributions import NegativeBinomial

        decoder_input = [z, y]
        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        if batch_index is not None:
            decoder_input.append(batch_index)

        h = self.decoder(*decoder_input)
        px_scale = self.px_decoder(h)
        px_rate = library * px_scale
        return {MODULE_KEYS.PX_KEY: NegativeBinomial(px_rate, logits=self.px_r)}

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution],
        generative_outputs: dict[str, Distribution],
        kl_weight: float = 1.0,
    ):
        """Loss computation."""
        from torch.distributions import Normal
        from torch.distributions import kl_divergence as kl

        from scvi.module.base import LossOutput

        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        qz = inference_outputs[MODULE_KEYS.QZ_KEY]
        px = generative_outputs[MODULE_KEYS.PX_KEY]

        mean = torch.zeros_like(qz.loc)
        scale = torch.ones_like(qz.scale)

        kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)

        reconst_loss = -px.log_prob(x).sum(-1)
        scaling_factor = self.ct_weight[y.long()[:, 0]]
        loss = torch.mean(scaling_factor * (reconst_loss + kl_weight * kl_divergence_z))

        return LossOutput(loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence_z)

    @torch.inference_mode()
    def sample(
        self,
        tensors: dict[str, torch.Tensor],
        n_samples: int = 1,
    ) -> torch.Tensor:
        r"""Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Tensors dict
        n_samples
            Number of required samples for each cell

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        inference_kwargs = {"n_samples": n_samples}
        generative_outputs = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )[1]

        dist = generative_outputs[MODULE_KEYS.PX_KEY]
        if n_samples > 1:
            exprs = dist.sample().permute([1, 2, 0])  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample()

        return exprs.cpu()
