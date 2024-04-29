from __future__ import annotations

from typing import Any

import flax.linen as nn
import jax
import jax.numpy as jnp
import numpy as np
import numpyro.distributions as dist

from scvi import REGISTRY_KEYS
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial
from scvi.external.mrvi._components import (
    AttentionBlock,
    ConditionalNormalization,
    Dense,
    NormalDistOutputNN,
)
from scvi.module.base import JaxBaseModuleClass, LossOutput, flax_configure

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

# Lower stddev leads to better initial loss values
_normal_initializer = jax.nn.initializers.normal(stddev=0.1)


class _DecoderZXAttention(nn.Module):
    n_in: int
    n_out: int
    n_batch: int
    n_latent_sample: int = 16
    h_activation: callable = nn.softmax
    n_channels: int = 4
    n_heads: int = 2
    dropout_rate: float = 0.1
    stop_gradients: bool = False
    stop_gradients_mlp: bool = False
    training: bool | None = None
    n_hidden: int = 32
    n_layers: int = 1
    training: bool | None = None
    low_dim_batch: bool = True
    activation: callable = nn.gelu

    @nn.compact
    def __call__(
        self,
        z: np.ndarray | jnp.ndarray,
        batch_covariate: np.ndarray | jnp.ndarray,
        size_factor: np.ndarray | jnp.ndarray,
        continuous_covariates: np.ndarray | jnp.ndarray | None,
        training: bool | None = None,
    ) -> NegativeBinomial:
        has_mc_samples = z.ndim == 3
        z_stop = z if not self.stop_gradients else jax.lax.stop_gradient(z)
        z_ = nn.LayerNorm(name="u_ln")(z_stop)

        batch_covariate = batch_covariate.astype(int).flatten()

        if self.n_batch >= 2:
            batch_embed = nn.Embed(
                self.n_batch, self.n_latent_sample, embedding_init=_normal_initializer
            )(batch_covariate)  # (batch, n_latent_sample)
            batch_embed = nn.LayerNorm(name="batch_embed_ln")(batch_embed)
            if has_mc_samples:
                batch_embed = jnp.tile(batch_embed, (z_.shape[0], 1, 1))

            res_dim = self.n_in if self.low_dim_batch else self.n_out

            query_embed = z_
            kv_embed = batch_embed
            residual = AttentionBlock(
                query_dim=self.n_in,
                out_dim=res_dim,
                outerprod_dim=self.n_latent_sample,
                n_channels=self.n_channels,
                n_heads=self.n_heads,
                dropout_rate=self.dropout_rate,
                n_hidden_mlp=self.n_hidden,
                n_layers_mlp=self.n_layers,
                stop_gradients_mlp=self.stop_gradients_mlp,
                training=training,
                activation=self.activation,
            )(query_embed=query_embed, kv_embed=kv_embed)

            if self.low_dim_batch:
                mu = nn.Dense(self.n_out)(z + residual)
            else:
                mu = nn.Dense(self.n_out)(z) + residual
        else:
            mu = nn.Dense(self.n_out)(z_)
        mu = self.h_activation(mu)
        return NegativeBinomial(
            mean=mu * size_factor,
            inverse_dispersion=jnp.exp(self.param("px_r", jax.random.normal, (self.n_out,))),
        )


class EncoderUZ(nn.Module):
    n_latent: int
    n_sample: int
    n_latent_u: int | None = None
    n_latent_sample: int = 16
    n_channels: int = 4
    n_heads: int = 2
    dropout_rate: float = 0.0
    stop_gradients: bool = False
    stop_gradients_mlp: bool = False
    use_map: bool = True
    n_hidden: int = 32
    n_layers: int = 1
    training: bool | None = None
    activation: callable = nn.gelu

    @nn.compact
    def __call__(
        self,
        u: np.ndarray | jnp.ndarray,
        sample_covariate: np.ndarray | jnp.ndarray,
        training: bool | None = None,
    ):
        training = nn.merge_param("training", self.training, training)
        sample_covariate = sample_covariate.astype(int).flatten()
        self.n_latent_u if self.n_latent_u is not None else self.n_latent  # noqa: B018
        has_mc_samples = u.ndim == 3
        u_stop = u if not self.stop_gradients else jax.lax.stop_gradient(u)
        u_ = nn.LayerNorm(name="u_ln")(u_stop)

        sample_embed = nn.Embed(
            self.n_sample, self.n_latent_sample, embedding_init=_normal_initializer
        )(sample_covariate)  # (batch, n_latent_sample)
        sample_embed = nn.LayerNorm(name="sample_embed_ln")(sample_embed)
        if has_mc_samples:
            sample_embed = jnp.tile(sample_embed, (u_.shape[0], 1, 1))

        n_outs = 1 if self.use_map else 2
        residual = AttentionBlock(
            query_dim=self.n_latent,
            out_dim=n_outs * self.n_latent,
            outerprod_dim=self.n_latent_sample,
            n_channels=self.n_channels,
            n_heads=self.n_heads,
            dropout_rate=self.dropout_rate,
            stop_gradients_mlp=self.stop_gradients_mlp,
            n_hidden_mlp=self.n_hidden,
            n_layers_mlp=self.n_layers,
            training=training,
            activation=self.activation,
        )(query_embed=u_, kv_embed=sample_embed)

        if self.n_latent_u is not None:
            z_base = nn.Dense(self.n_latent)(u_stop)
            return z_base, residual
        else:
            return u, residual


class _EncoderXU(nn.Module):
    n_latent: int
    n_sample: int
    n_hidden: int
    n_layers: int = 1
    activation: callable = nn.gelu
    training: bool | None = None

    @nn.compact
    def __call__(
        self,
        x: np.ndarray | jnp.ndarray,
        sample_covariate: np.ndarray | jnp.ndarray,
        training: bool | None = None,
    ) -> dist.Normal:
        training = nn.merge_param("training", self.training, training)
        x_feat = jnp.log1p(x)
        for _ in range(2):
            x_feat = Dense(self.n_hidden)(x_feat)
            x_feat = ConditionalNormalization(self.n_hidden, self.n_sample)(
                x_feat, sample_covariate, training=training
            )
            x_feat = self.activation(x_feat)
        sample_effect = nn.Embed(self.n_sample, self.n_hidden, embedding_init=_normal_initializer)(
            sample_covariate.squeeze(-1).astype(int)
        )
        inputs = x_feat + sample_effect
        return NormalDistOutputNN(self.n_latent, self.n_hidden, self.n_layers)(
            inputs, training=training
        )


@flax_configure
class MrVAE(JaxBaseModuleClass):
    """Flax module for the Multi-resolution Variational Inference (MrVI) model."""

    n_input: int
    n_sample: int
    n_batch: int
    n_labels: int
    n_continuous_cov: int
    n_latent: int = 30
    n_latent_u: int = 10
    encoder_n_hidden: int = 128
    encoder_n_layers: int = 2
    z_u_prior: bool = True
    z_u_prior_scale: float = 0.0
    u_prior_scale: float = 0.0
    u_prior_mixture: bool = True
    u_prior_mixture_k: int = 20
    learn_z_u_prior_scale: bool = False
    laplace_scale: float = None
    scale_observations: bool = False
    px_kwargs: dict | None = None
    qz_kwargs: dict | None = None
    qu_kwargs: dict | None = None
    training: bool = True
    n_obs_per_sample: jnp.ndarray | None = None

    def setup(self):
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

        # Generative model
        px_cls = _DecoderZXAttention
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
        self.qu = _EncoderXU(
            n_latent=self.n_latent if is_isomorphic_uz else n_latent_u,
            n_sample=self.n_sample,
            n_hidden=self.encoder_n_hidden,
            n_layers=self.encoder_n_layers,
            **qu_kwargs,
        )

        if self.learn_z_u_prior_scale:
            self.pz_scale = self.param("pz_scale", nn.initializers.zeros, (self.n_latent,))
        else:
            self.pz_scale = self.z_u_prior_scale

        if self.u_prior_mixture:
            if self.n_labels > 1:
                u_prior_mixture_k = self.n_labels
            else:
                u_prior_mixture_k = self.u_prior_mixture_k
            u_dim = self.n_latent_u if self.n_latent_u is not None else self.n_latent
            self.u_prior_logits = self.param(
                "u_prior_logits", nn.initializers.zeros, (u_prior_mixture_k,)
            )
            self.u_prior_means = self.param(
                "u_prior_means", jax.random.normal, (u_dim, u_prior_mixture_k)
            )
            self.u_prior_scales = self.param(
                "u_prior_scales", nn.initializers.zeros, (u_dim, u_prior_mixture_k)
            )

    @property
    def required_rngs(self):
        return ("params", "u", "dropout", "eps")

    def _get_inference_input(self, tensors: dict[str, np.ndarray | jnp.ndarray]) -> dict[str, Any]:
        x = tensors[REGISTRY_KEYS.X_KEY]
        sample_index = tensors[REGISTRY_KEYS.SAMPLE_KEY]
        return {"x": x, "sample_index": sample_index}

    def inference(self, x, sample_index, mc_samples=None, cf_sample=None, use_mean=False):
        """Latent variable inference."""
        qu = self.qu(x, sample_index, training=self.training)
        if use_mean:
            u = qu.mean
        else:
            u_rng = self.make_rng("u")
            sample_shape = (mc_samples,) if mc_samples is not None else ()
            u = qu.rsample(u_rng, sample_shape=sample_shape)

        sample_index_cf = sample_index if cf_sample is None else cf_sample

        z_base, eps = self.qz(u, sample_index_cf, training=self.training)
        qeps_ = eps

        qeps = None
        if qeps_.shape[-1] == 2 * self.n_latent:
            loc_, scale_ = qeps_[..., : self.n_latent], qeps_[..., self.n_latent :]
            qeps = dist.Normal(loc_, nn.softplus(scale_) + 1e-3)
            eps = qeps.mean if use_mean else qeps.rsample(self.make_rng("eps"))
        As = None
        z = z_base + eps
        library = jnp.log(x.sum(1, keepdims=True))

        return {
            "qu": qu,
            "qeps": qeps,
            "eps": eps,
            "As": As,
            "u": u,
            "z": z,
            "z_base": z_base,
            "library": library,
        }

    def _get_generative_input(
        self,
        tensors: dict[str, np.ndarray | jnp.ndarray],
        inference_outputs: dict[str, Any],
    ) -> dict[str, Any]:
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        label_index = tensors[REGISTRY_KEYS.LABELS_KEY]
        continuous_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None)
        return {
            "z": z,
            "library": library,
            "batch_index": batch_index,
            "label_index": label_index,
            "continuous_covs": continuous_covs,
        }

    def generative(self, z, library, batch_index, label_index, continuous_covs):
        """Generative model."""
        library_exp = jnp.exp(library)
        px = self.px(
            z,
            batch_index,
            size_factor=library_exp,
            continuous_covariates=continuous_covs,
            training=self.training,
        )
        h = px.mean / library_exp

        if self.u_prior_mixture:
            offset = (
                10.0 * jax.nn.one_hot(label_index, self.n_labels) if self.n_labels >= 2 else 0.0
            )
            cats = dist.Categorical(logits=self.u_prior_logits + offset)
            normal_dists = dist.Normal(self.u_prior_means, jnp.exp(self.u_prior_scales))
            pu = dist.MixtureSameFamily(cats, normal_dists)
        else:
            pu = dist.Normal(0, jnp.exp(self.u_prior_scale))
        return {"px": px, "pu": pu, "h": h}

    def loss(
        self,
        tensors: dict[str, np.ndarray | jnp.ndarray],
        inference_outputs: dict[str, Any],
        generative_outputs: dict[str, Any],
        kl_weight: float = 1.0,
    ) -> jnp.ndarray:
        """Compute the loss function value."""
        reconstruction_loss = (
            -generative_outputs["px"].log_prob(tensors[REGISTRY_KEYS.X_KEY]).sum(-1)
        )

        if self.u_prior_mixture:
            kl_u = inference_outputs["qu"].log_prob(inference_outputs["u"]) - generative_outputs[
                "pu"
            ].log_prob(inference_outputs["u"])
            kl_u = kl_u.sum(-1)
        else:
            kl_u = dist.kl_divergence(inference_outputs["qu"], generative_outputs["pu"]).sum(-1)
        inference_outputs["qeps"]

        kl_z = 0.0
        eps = inference_outputs["z"] - inference_outputs["z_base"]
        if self.z_u_prior:
            peps = dist.Normal(0, jnp.exp(self.pz_scale))
            kl_z = -peps.log_prob(eps).sum(-1)

        weighted_kl_local = kl_weight * (kl_u + kl_z)
        loss = reconstruction_loss + weighted_kl_local

        if self.laplace_scale is not None:
            As = inference_outputs["As"]
            n_obs = As.shape[0]
            As = As.reshape((n_obs, -1))
            p_As = dist.Laplace(0, self.laplace_scale).log_prob(As).sum(-1)

            n_obs_total = self.n_obs_per_sample.sum()
            As_pen = -p_As / n_obs_total
            As_pen = As_pen.sum()
            loss = loss + (kl_weight * As_pen)

        if self.scale_observations:
            sample_index = tensors[REGISTRY_KEYS.SAMPLE_KEY].flatten().astype(int)
            prefactors = self.n_obs_per_sample[sample_index]
            loss = loss / prefactors

        loss = jnp.mean(loss)

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            kl_local=(kl_u + kl_z),
        )

    def compute_h_from_x_eps(
        self,
        x,
        sample_index,
        batch_index,
        extra_eps,
        cf_sample=None,
        continuous_covs=None,
        mc_samples=10,
    ):
        """Compute normalized gene expression from observations using predefined eps"""
        library = 7.0 * jnp.ones_like(
            sample_index
        )  # placeholder, has no effect on the value of h.
        inference_outputs = self.inference(
            x, sample_index, mc_samples=mc_samples, cf_sample=cf_sample, use_mean=False
        )
        generative_inputs = {
            "z": inference_outputs["z_base"] + extra_eps,
            "library": library,
            "batch_index": batch_index,
            "continuous_covs": continuous_covs,
            "label_index": jnp.zeros([x.shape[0], 1]),
        }
        generative_outputs = self.generative(**generative_inputs)
        return generative_outputs["h"]
