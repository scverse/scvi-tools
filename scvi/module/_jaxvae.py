from typing import Dict, Optional

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial
from scvi.module.base import JaxBaseModuleClass, LossOutput, flax_configure


class Dense(nn.Dense):
    """Jax dense layer."""

    def __init__(self, *args, **kwargs):
        # scale set to reimplement pytorch init
        scale = 1 / 3
        kernel_init = variance_scaling(scale, "fan_in", "uniform")
        # bias init can't see input shape so don't include here
        kwargs.update({"kernel_init": kernel_init})
        super().__init__(*args, **kwargs)


class FlaxEncoder(nn.Module):
    """Encoder for Jax VAE."""

    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int
    training: Optional[bool] = None

    def setup(self):
        """Setup encoder."""
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_latent)
        self.dense4 = Dense(self.n_latent)

        self.batchnorm1 = nn.BatchNorm(momentum=0.99, epsilon=0.001)
        self.batchnorm2 = nn.BatchNorm(momentum=0.99, epsilon=0.001)
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

    def __call__(self, x: jnp.ndarray, training: Optional[bool] = None):
        """Forward pass."""
        training = nn.merge_param("training", self.training, training)
        is_eval = not training

        x_ = jnp.log1p(x)

        h = self.dense1(x_)
        h = self.batchnorm1(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense2(h)
        h = self.batchnorm2(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout2(h, deterministic=is_eval)

        mean = self.dense3(h)
        log_var = self.dense4(h)

        return mean, jnp.exp(log_var)


class FlaxDecoder(nn.Module):
    """Decoder for Jax VAE."""

    n_input: int
    dropout_rate: float
    n_hidden: int
    training: Optional[bool] = None

    def setup(self):
        """Setup decoder."""
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_hidden)
        self.dense4 = Dense(self.n_hidden)
        self.dense5 = Dense(self.n_input)

        self.batchnorm1 = nn.BatchNorm(momentum=0.99, epsilon=0.001)
        self.batchnorm2 = nn.BatchNorm(momentum=0.99, epsilon=0.001)
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

        self.disp = self.param(
            "disp", lambda rng, shape: jax.random.normal(rng, shape), (self.n_input, 1)
        )

    def __call__(
        self, z: jnp.ndarray, batch: jnp.ndarray, training: Optional[bool] = None
    ):
        """Forward pass."""
        # TODO(adamgayoso): Test this
        training = nn.merge_param("training", self.training, training)
        is_eval = not training

        h = self.dense1(z)
        h += self.dense2(batch)

        h = self.batchnorm1(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense3(h)
        # skip connection
        h += self.dense4(batch)
        h = self.batchnorm2(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout2(h, deterministic=is_eval)
        h = self.dense5(h)
        return h, self.disp.ravel()


@flax_configure
class JaxVAE(JaxBaseModuleClass):
    """Variational autoencoder model."""

    n_input: int
    n_batch: int
    n_hidden: Tunable[int] = 128
    n_latent: Tunable[int] = 30
    dropout_rate: Tunable[float] = 0.0
    n_layers: Tunable[int] = 1
    gene_likelihood: Tunable[str] = "nb"
    eps: Tunable[float] = 1e-8
    training: bool = True

    def setup(self):
        """Setup model."""
        self.encoder = FlaxEncoder(
            n_input=self.n_input,
            n_latent=self.n_latent,
            n_hidden=self.n_hidden,
            dropout_rate=self.dropout_rate,
        )

        self.decoder = FlaxDecoder(
            n_input=self.n_input,
            dropout_rate=0.0,
            n_hidden=self.n_hidden,
        )

    @property
    def required_rngs(self):  # noqa: D102
        return ("params", "dropout", "z")

    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray]):
        """Get input for inference."""
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = dict(x=x)
        return input_dict

    def inference(self, x: jnp.ndarray, n_samples: int = 1) -> dict:
        """Run inference model."""
        mean, var = self.encoder(x, training=self.training)
        stddev = jnp.sqrt(var) + self.eps

        qz = dist.Normal(mean, stddev)
        z_rng = self.make_rng("z")
        sample_shape = () if n_samples == 1 else (n_samples,)
        z = qz.rsample(z_rng, sample_shape=sample_shape)

        return dict(qz=qz, z=z)

    def _get_generative_input(
        self,
        tensors: Dict[str, jnp.ndarray],
        inference_outputs: Dict[str, jnp.ndarray],
    ):
        """Get input for generative model."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        z = inference_outputs["z"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = dict(
            x=x,
            z=z,
            batch_index=batch_index,
        )
        return input_dict

    def generative(self, x, z, batch_index) -> dict:
        """Run generative model."""
        # one hot adds an extra dimension
        batch = jax.nn.one_hot(batch_index, self.n_batch).squeeze(-2)
        rho_unnorm, disp = self.decoder(z, batch, training=self.training)
        disp_ = jnp.exp(disp)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho

        if self.gene_likelihood == "nb":
            disp_ = jnp.exp(disp)
            px = NegativeBinomial(mean=mu, inverse_dispersion=disp_)
        else:
            px = dist.Poisson(mu)

        return dict(px=px, rho=rho)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Compute loss."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        px = generative_outputs["px"]
        qz = inference_outputs["qz"]
        reconst_loss = -px.log_prob(x).sum(-1)
        kl_divergence_z = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)

        kl_local_for_warmup = kl_divergence_z
        weighted_kl_local = kl_weight * kl_local_for_warmup

        loss = jnp.mean(reconst_loss + weighted_kl_local)

        kl_local = kl_divergence_z
        return LossOutput(
            loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_local
        )
