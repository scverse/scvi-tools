from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import torch
import torch.distributions as dist
import torch.nn as nn
import torch.nn.init as init

from scvi import REGISTRY_KEYS, settings
from scvi.distributions import NegativeBinomial
from scvi.external.mrvi_torch._components import (
    AttentionBlock,
    ConditionalNormalization,
    NormalDistOutputNN,
)
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

if TYPE_CHECKING:
    from collections.abc import Callable

DEFAULT_PX_KWARGS = {
    "n_hidden": 32,
    "stop_gradients": False,
    "stop_gradients_mlp": True,
    "dropout_rate": 0.03,
}
DEFAULT_QZ_ATTENTION_KWARGS = {
    "use_map": True,
    "stop_gradients": False,
    "stop_gradients_mlp": True,
    "dropout_rate": 0.03,
}
DEFAULT_QU_KWARGS = {}


class DecoderZXAttention(nn.Module):
    """Attention-based decoder.

    Parameters
    ----------
    n_in
        Number of input features.
    n_out
        Number of output features.
    n_batch
        Number of batches.
    n_latent_sample
        Number of latent samples.
    h_activation
        Activation function for the output layer.
    n_channels
        Number of channels in the attention block.
    n_heads
        Number of heads in the attention block.
    dropout_rate
        Dropout rate.
    stop_gradients
        Whether to stop gradients to ``z``.
    stop_gradients_mlp
        Whether to stop gradients to the MLP in the attention block.
    n_hidden
        Number of hidden units in the MLP.
    n_layers
        Number of layers in the MLP.
    low_dim_batch
        Whether to use low-dimensional batch embeddings.
    activation
        Activation function for the MLP.
    """

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_batch: int,
        n_latent_sample: int = 16,
        h_activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.softmax,
        n_channels: int = 4,
        n_heads: int = 2,
        dropout_rate: float = 0.1,
        stop_gradients: bool = False,
        stop_gradients_mlp: bool = False,
        n_hidden: int = 32,
        n_layers: int = 1,
        low_dim_batch: bool = True,
        activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.gelu,
    ):
        super().__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.n_batch = n_batch
        self.n_latent_sample = n_latent_sample
        self.h_activation = h_activation
        self.n_channels = n_channels
        self.n_heads = n_heads
        self.dropout_rate = dropout_rate
        self.stop_gradients = stop_gradients
        self.stop_gradients_mlp = stop_gradients_mlp
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.low_dim_batch = low_dim_batch
        self.activation = activation

        self.layer_norm = nn.LayerNorm(self.n_in)
        self.layer_norm_batch_embed = nn.LayerNorm(self.n_latent_sample)
        self.batch_embedding = nn.Embedding(self.n_batch, self.n_latent_sample)
        # Lower stddev leads to better initial loss values
        init.normal_(self.batch_embedding.weight, std=0.1)

        self.res_dim = self.n_in if self.low_dim_batch else self.n_out

        self.attention_block = AttentionBlock(
            query_dim=self.n_in,
            kv_dim=self.n_latent_sample,
            out_dim=self.res_dim,
            outerprod_dim=self.n_latent_sample,
            n_channels=self.n_channels,
            n_heads=self.n_heads,
            dropout_rate=self.dropout_rate,
            n_hidden_mlp=self.n_hidden,
            n_layers_mlp=self.n_layers,
            stop_gradients_mlp=self.stop_gradients_mlp,
            activation=self.activation,
        )

        self.fc = nn.Linear(self.n_in, self.n_out)

        self.px_r = nn.Parameter(
            torch.zeros(
                self.n_out,
            ),
            requires_grad=True,
        )
        init.normal_(self.px_r)

    def forward(
        self,
        z: torch.Tensor,
        batch_covariate: torch.Tensor,
        size_factor: torch.Tensor,
    ) -> NegativeBinomial:
        has_mc_samples = z.ndim == 3
        z_stop = z if not self.stop_gradients else z.detach()
        z_ = self.layer_norm(z_stop)

        batch_covariate = batch_covariate.to(torch.int64).flatten()

        if self.n_batch >= 2:
            batch_embed = self.batch_embedding(batch_covariate)
            batch_embed = self.layer_norm_batch_embed(batch_embed)
            if has_mc_samples:
                batch_embed = batch_embed.unsqueeze(0).expand(z_.shape[0], -1, -1)

            query_embed = z_
            kv_embed = batch_embed
            residual = self.attention_block(query_embed=query_embed, kv_embed=kv_embed)
            if self.low_dim_batch:
                mu = self.fc(z + residual)
            else:
                mu = self.fc(z) + residual
        else:
            mu = self.fc(z_)
        if self.h_activation in [nn.functional.softmax, torch.softmax]:
            mu = self.h_activation(mu, dim=-1)
        else:
            mu = self.h_activation(mu)
        return NegativeBinomial(
            mu=mu * size_factor,
            theta=torch.exp(self.px_r),
        )


class EncoderUZ(nn.Module):
    """Attention-based encoder from ``u`` to ``z``.

    Parameters
    ----------
    n_latent
        Number of latent variables.
    n_sample
        Number of samples.
    n_latent_u
        Number of latent variables for ``u``.
    n_latent_sample
        Number of latent samples.
    n_channels
        Number of channels in the attention block.
    n_heads
        Number of heads in the attention block.
    dropout_rate
        Dropout rate.
    stop_gradients
        Whether to stop gradients to ``u``.
    stop_gradients_mlp
        Whether to stop gradients to the MLP in the attention block.
    use_map
        Whether to use the MAP estimate to approximate the posterior of ``z`` given ``u``
    n_hidden
        Number of hidden units in the MLP.
    n_layers
        Number of layers in the MLP.
    activation
        Activation function for the MLP.
    """

    def __init__(
        self,
        n_latent: int,
        n_sample: int,
        n_latent_u: int | None = None,
        n_latent_sample: int = 16,
        n_channels: int = 4,
        n_heads: int = 2,
        dropout_rate: float = 0.0,
        stop_gradients: bool = False,
        stop_gradients_mlp: bool = False,
        use_map: bool = True,
        n_hidden: int = 32,
        n_layers: int = 1,
        activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.gelu,
    ):
        super().__init__()
        self.n_latent = n_latent
        self.n_sample = n_sample
        self.n_latent_u = n_latent_u if n_latent_u is not None else n_latent
        self.n_latent_sample = n_latent_sample
        self.n_channels = n_channels
        self.n_heads = n_heads
        self.dropout_rate = dropout_rate
        self.stop_gradients = stop_gradients
        self.stop_gradients_mlp = stop_gradients_mlp
        self.use_map = use_map
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.activation = activation

        self.layer_norm = nn.LayerNorm(self.n_latent_u)
        self.embedding = nn.Embedding(self.n_sample, self.n_latent_sample)
        # Initialize with same standard deviation as JAX version
        init.normal_(self.embedding.weight, std=0.1)
        self.layer_norm_embed = nn.LayerNorm(self.n_latent_sample)

        n_outs = 1 if self.use_map else 2
        self.attention_block = AttentionBlock(
            query_dim=self.n_latent_u,
            kv_dim=self.n_latent_sample,
            out_dim=n_outs * self.n_latent,
            outerprod_dim=self.n_latent_sample,
            n_channels=self.n_channels,
            n_heads=self.n_heads,
            dropout_rate=self.dropout_rate,
            stop_gradients_mlp=self.stop_gradients_mlp,
            n_hidden_mlp=self.n_hidden,
            n_layers_mlp=self.n_layers,
            activation=self.activation,
        )

        if self.n_latent_u is not None:
            self.fc = nn.Linear(self.n_latent_u, self.n_latent)

    def forward(
        self,
        u: torch.Tensor,
        sample_covariate: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        sample_covariate = sample_covariate.to(torch.int64).flatten()
        has_mc_samples = u.ndim == 3
        u_stop = u if not self.stop_gradients else u.detach()
        u_ = self.layer_norm(u_stop)

        sample_embed = self.layer_norm_embed(self.embedding(sample_covariate))
        # nn.init.normal_(self.sample_embed.weight, mean=0.0, std=0.1)

        if has_mc_samples:
            sample_embed = sample_embed.repeat(u_.shape[0], 1, 1)

        residual = self.attention_block(query_embed=u_, kv_embed=sample_embed)

        if self.n_latent_u is not None:
            z_base = self.fc(u_stop)
            return z_base, residual
        else:
            return u, residual


class EncoderXU(nn.Module):
    """Encoder from ``x`` to ``u``.

    Parameters
    ----------
    n_input
        Number of input variables.
    n_latent
        Number of latent variables.
    n_sample
        Number of samples.
    n_hidden
        Number of hidden units in the MLP.
    n_layers
        Number of layers in the MLP.
    training
        Whether the model is in training mode.
    """

    def __init__(
        self,
        n_input: int,
        n_latent: int,
        n_sample: int,
        n_hidden: int = 128,
        n_layers: int = 1,
        activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.gelu,
        training: bool | None = None,
    ):
        super().__init__()
        self.n_input = n_input
        self.n_latent = n_latent
        self.n_sample = n_sample
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.activation = activation
        self.training = training

        self.fc1 = nn.Linear(self.n_input, self.n_hidden)
        self.conditional_norm1 = ConditionalNormalization(self.n_hidden, self.n_sample)
        self.conditional_norm1.train()
        self.fc2 = nn.Linear(self.n_hidden, self.n_hidden)
        self.conditional_norm2 = ConditionalNormalization(self.n_hidden, self.n_sample)
        self.conditional_norm2.train()

        self.sample_embed = nn.Embedding(self.n_sample, self.n_hidden)
        # Initialize with same standard deviation as JAX version
        init.normal_(self.sample_embed.weight, std=0.1)

        # Create the output network in __init__ instead of forward
        self.output_nn = NormalDistOutputNN(
            self.n_hidden, self.n_latent, self.n_hidden, self.n_layers
        )

    def forward(
        self, x: torch.Tensor, sample_covariate: torch.Tensor, training: bool | None = None
    ) -> dist.Normal:
        training = training if training is not None else self.training
        x_feat = torch.log1p(x)
        x_feat = self.fc1(x_feat)
        x_feat = self.conditional_norm1(x_feat, sample_covariate, training=training)
        x_feat = self.activation(x_feat)
        x_feat = self.fc2(x_feat)
        x_feat = self.conditional_norm2(x_feat, sample_covariate, training=training)
        x_feat = self.activation(x_feat)
        sample_effect = self.sample_embed(sample_covariate.squeeze(-1).to(torch.int64))
        inputs = x_feat + sample_effect
        return self.output_nn(inputs)


class TorchMRVAE(BaseModuleClass):
    """Multi-resolution Variational Inference (MrVI) module.

    Parameters
    ----------
    n_input
        Number of input features.
    n_sample
        Number of samples.
    n_batch
        Number of batches.
    n_labels
        Number of labels.
    n_latent
        Number of latent variables for ``z``.
    n_latent_u
        Number of latent variables for ``u``.
    encoder_n_hidden
        Number of hidden units in the encoder.
    encoder_n_layers
        Number of layers in the encoder.
    z_u_prior
        Whether to place a Gaussian prior on ``z`` given ``u``.
    z_u_prior_scale
        Natural log of the scale parameter of the Gaussian prior placed on ``z`` given ``u``. Only
        applies of ``learn_z_u_prior_scale`` is ``False``.
    u_prior_scale
        Natural log of the scale parameter of the Gaussian prior placed on ``u``. If
        ``u_prior_mixture`` is ``True``, this scale applies to each mixture component distribution.
    u_prior_mixture
        Whether to use a mixture of Gaussians prior for ``u``.
    u_prior_mixture_k
        Number of mixture components to use for the mixture of Gaussians prior on ``u``.
    learn_z_u_prior_scale
        Whether to learn the scale parameter of the prior distribution of ``z`` given ``u``.
    scale_observations
        Whether to scale the loss associated with each observation by the total number of
        observations linked to the associated sample.
    px_kwargs
        Keyword arguments for the generative model.
    qz_kwargs
        Keyword arguments for the inference model from ``u`` to ``z``.
    qu_kwargs
        Keyword arguments for the inference model from ``x`` to ``u``.
    training
        Whether the model is in training mode.
    n_obs_per_sample
        Number of observations per sample.
    """

    def __init__(
        self,
        n_input: int,
        n_sample: int,
        n_batch: int,
        n_labels: int,
        n_latent: int = 30,
        n_latent_u: int = 10,
        encoder_n_hidden: int = 128,
        encoder_n_layers: int = 2,
        z_u_prior: bool = True,
        z_u_prior_scale: torch.Tensor = torch.zeros(1, dtype=torch.float32),
        u_prior_scale: float = 0.0,
        u_prior_mixture: bool = True,
        u_prior_mixture_k: int = 20,
        learn_z_u_prior_scale: bool = False,
        scale_observations: bool = False,
        px_kwargs: dict | None = None,
        qz_kwargs: dict | None = None,
        qu_kwargs: dict | None = None,
        n_obs_per_sample: torch.Tensor | None = None,
    ):
        super().__init__()
        self.n_input = n_input
        self.n_sample = n_sample
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_latent = n_latent
        self.n_latent_u = n_latent_u
        self.encoder_n_hidden = encoder_n_hidden
        self.encoder_n_layers = encoder_n_layers
        self.z_u_prior = z_u_prior
        self.z_u_prior_scale = z_u_prior_scale
        self.u_prior_scale = u_prior_scale
        self.u_prior_mixture = u_prior_mixture
        self.u_prior_mixture_k = u_prior_mixture_k
        self.learn_z_u_prior_scale = learn_z_u_prior_scale
        self.scale_observations = scale_observations
        self.px_kwargs = px_kwargs
        self.qz_kwargs = qz_kwargs
        self.qu_kwargs = qu_kwargs
        self.n_obs_per_sample = n_obs_per_sample

        px_kwargs = DEFAULT_PX_KWARGS.copy()
        if self.px_kwargs is not None:
            px_kwargs.update(self.px_kwargs)

        qz_kwargs = DEFAULT_QZ_ATTENTION_KWARGS.copy()
        if self.qz_kwargs is not None:
            qz_kwargs.update(self.qz_kwargs)

        qu_kwargs = DEFAULT_QU_KWARGS.copy()
        if self.qu_kwargs is not None:
            qu_kwargs.update(self.qu_kwargs)

        is_isomorphic_uz = self.n_latent == self.n_latent_u
        n_latent_u = None if is_isomorphic_uz else self.n_latent_u

        if self.n_latent < self.n_latent_u:
            warnings.warn(
                "The number of latent variables for `z` is set to less than the number of latent "
                "variables for `u`.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        # Generative model
        px_cls = DecoderZXAttention
        self.px = px_cls(
            self.n_latent,
            self.n_input,
            self.n_batch,
            **px_kwargs,
        )

        qz_cls = EncoderUZ
        self.qz = qz_cls(
            self.n_latent,
            self.n_sample,
            n_latent_u=n_latent_u,
            **qz_kwargs,
        )

        # Inference model
        self.qu = EncoderXU(
            n_input=self.n_input,
            n_latent=self.n_latent if is_isomorphic_uz else n_latent_u,
            n_sample=self.n_sample,
            n_hidden=self.encoder_n_hidden,
            n_layers=self.encoder_n_layers,
            **qu_kwargs,
        )
        self.backend = "torch"

        if self.learn_z_u_prior_scale:
            self.pz_scale = nn.Parameter(torch.zeros(self.n_latent))
        else:
            self.pz_scale = self.z_u_prior_scale

        if self.u_prior_mixture:
            if self.n_labels > 1:
                u_prior_mixture_k = self.n_labels
            else:
                u_prior_mixture_k = self.u_prior_mixture_k

            u_dim = self.n_latent_u if self.n_latent_u is not None else self.n_latent
            self.u_prior_logits = nn.Parameter(torch.zeros(u_prior_mixture_k))
            self.u_prior_means = nn.Parameter(torch.randn(u_prior_mixture_k, u_dim))
            self.u_prior_scales = nn.Parameter(torch.zeros(u_prior_mixture_k, u_dim))

    def _get_inference_input(self, tensors: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
        x = tensors[REGISTRY_KEYS.X_KEY]
        sample_index = tensors[REGISTRY_KEYS.SAMPLE_KEY]
        return {"x": x, "sample_index": sample_index}

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        sample_index: torch.Tensor,
        mc_samples: int | None = None,
        cf_sample: torch.Tensor | None = None,
        use_mean: bool = False,
    ) -> dict[str, torch.Tensor]:
        if type(x).__name__ != "Tensor" and x is not None:
            x = torch.Tensor(x)
        if type(sample_index).__name__ != "Tensor" and sample_index is not None:
            sample_index = torch.Tensor(sample_index)
        qu = self.qu(x, sample_index, training=self.training)
        if use_mean:
            u = qu.mean
        else:
            sample_shape = (mc_samples,) if mc_samples is not None else ()
            # Use functional random sampling instead of rsample to avoid vmap issues
            noise_shape = sample_shape + qu.mean.shape
            noise = torch.normal(
                mean=0.0, std=1.0, size=noise_shape, device=qu.mean.device, dtype=qu.mean.dtype
            )
            u = qu.mean + noise * qu.scale

        sample_index_cf = sample_index if cf_sample is None else cf_sample

        z_base, eps = self.qz(u, sample_index_cf)
        qeps_ = eps

        qeps = None
        if qeps_.shape[-1] == 2 * self.n_latent:
            loc_, scale_ = qeps_[..., : self.n_latent], qeps_[..., self.n_latent :]
            qeps = dist.Normal(loc_, nn.functional.softplus(scale_) + 1e-3)
            if use_mean:
                eps = qeps.mean
            else:
                # Use functional random sampling instead of rsample to avoid vmap issues
                noise = torch.normal(
                    mean=0.0,
                    std=1.0,
                    size=qeps.mean.shape,
                    device=qeps.mean.device,
                    dtype=qeps.mean.dtype,
                )
                eps = qeps.mean + noise * qeps.scale
        z = z_base + eps
        library = torch.log(x.sum(1, keepdims=True))

        return {
            "qu": qu,
            "qeps": qeps,
            "eps": eps,
            "u": u,
            "z": z,
            "z_base": z_base,
            "library": library,
        }

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
    ) -> dict[str, torch.Tensor]:
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        label_index = tensors[REGISTRY_KEYS.LABELS_KEY]
        return {
            "z": z,
            "library": library,
            "batch_index": batch_index,
            "label_index": label_index,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor,
        label_index: torch.Tensor,
    ) -> dict[str, torch.Tensor | dist.Distribution]:
        """Generative model."""
        library_exp = torch.exp(library)
        px = self.px(z, batch_index, size_factor=library_exp)
        h = px.mean / library_exp

        if self.u_prior_mixture:
            offset = (
                10.0 * nn.functional.one_hot(label_index, self.n_labels)
                if self.n_labels >= 2
                else 0.0
            )
            cats = dist.Categorical(logits=self.u_prior_logits + offset)

            # TODO: do I need to use the independent as in below and why? Does jax version do this?
            normal_dists = dist.Independent(
                dist.Normal(self.u_prior_means, torch.exp(self.u_prior_scales)), 1
            )

            pu = dist.MixtureSameFamily(cats, normal_dists)
        else:
            pu = dist.Normal(0, torch.exp(self.u_prior_scale))
        return {"px": px, "pu": pu, "h": h}

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        kl_weight: float = 1.0,
    ) -> LossOutput:
        """Compute the loss function value."""
        reconstruction_loss = (
            -generative_outputs["px"].log_prob(tensors[REGISTRY_KEYS.X_KEY]).sum(-1)
        )

        if self.u_prior_mixture:
            kl_u = inference_outputs["qu"].log_prob(inference_outputs["u"]).sum(
                -1
            ) - generative_outputs["pu"].log_prob(inference_outputs["u"])
        else:
            kl_u = dist.kl_divergence(inference_outputs["qu"], generative_outputs["pu"]).sum(-1)

        kl_z = 0.0
        eps = inference_outputs["z"] - inference_outputs["z_base"]
        device = self.pz_scale.device
        eps = eps.to(device)
        if self.z_u_prior:
            peps = dist.Normal(0, torch.exp(self.pz_scale))
            kl_z = -peps.log_prob(eps).sum(-1)
            kl_z = kl_z.to(kl_u.device)

        kl_weight = (
            kl_weight.to(kl_u.device)
            if isinstance(kl_weight, torch.Tensor)
            else torch.tensor(kl_weight, device=kl_u.device)
        )
        weighted_kl_local = kl_weight * (kl_u + kl_z)
        loss = reconstruction_loss + weighted_kl_local

        if self.scale_observations:
            sample_index = tensors[REGISTRY_KEYS.SAMPLE_KEY].flatten().to(torch.int64)
            prefactors = self.n_obs_per_sample.to(self.device)[sample_index]
            loss = loss / prefactors

        loss = torch.mean(loss)

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            kl_local=(kl_u + kl_z),
        )

    def compute_h_from_x_eps(
        self,
        x: torch.Tensor,
        sample_index: torch.Tensor,
        batch_index: torch.Tensor,
        extra_eps: float,
        cf_sample: torch.Tensor | None = None,
        mc_samples: int = 10,
    ):
        """Compute normalized gene expression from observations using predefined eps"""
        library = 7.0 * torch.ones_like(
            sample_index
        )  # placeholder, has no effect on the value of h.
        inference_outputs = self.inference(
            x, sample_index, mc_samples=mc_samples, cf_sample=cf_sample, use_mean=False
        )
        generative_inputs = {
            "z": inference_outputs["z_base"] + extra_eps,
            "library": library,
            "batch_index": batch_index,
            "label_index": torch.zeros([x.shape[0], 1]),
        }
        generative_outputs = self.generative(**generative_inputs)
        return generative_outputs["h"]
