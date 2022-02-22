from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn

from scvi import REGISTRY_KEYS


class FlaxEncoder(nn.Module):
    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int

    @nn.compact
    def __call__(self, inputs, is_training):
        inputs_ = jnp.log1p(inputs)

        # k1, k2 = hk.next_rng_keys(2)

        h = nn.Dense(self.n_hidden)(inputs_)
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        # h = nn.Dropout(dropout_rate)(h)
        h = nn.Dense(self.n_hidden)(h)
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        # h = nn.Dropout(dropout_rate)(h)

        mean = nn.Dense(self.n_latent)(h)
        log_concentration = nn.Dense(self.n_latent)(h)

        return mean, nn.softplus(log_concentration)


class FlaxDecoder(nn.Module):
    n_input: int
    dropout_rate: float
    n_hidden: int

    @nn.compact
    def __call__(self, inputs, is_training):
        # disp = nn.get_parameter("disp", (self.n_input,1), init=jnp.ones)
        disp = self.param("disp", lambda rng, shape: jnp.ones(shape), (self.n_input, 1))

        # k1, k2 = nn.next_rng_keys(2)

        # print(inputs.shape, self.n_hidden)
        h = nn.Dense(self.n_hidden)(inputs)
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        # h = nn.Dropout(dropout_rate)(h)
        # skip connection
        h = nn.Dense(self.n_hidden)(jnp.concatenate([h, inputs], axis=-1))
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        # h = nn.Dropout(dropout_rate)(h)
        h = nn.Dense(self.n_input)(h)
        return h, disp.ravel()


class VAEOutput(NamedTuple):
    mean: jnp.ndarray
    stddev: jnp.ndarray
    nb: dist.NegativeBinomialLogits


class JaxVAE(nn.Module):
    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    is_training: bool = False
    n_layers: int = 1

    @nn.compact
    def __call__(self, array_dict, z_rng) -> VAEOutput:

        x = array_dict[REGISTRY_KEYS.X_KEY]
        batch = array_dict[REGISTRY_KEYS.BATCH_KEY]

        # one hot adds an extra dimension
        batch = jax.nn.one_hot(batch, self.n_batch).squeeze(-2)
        mean, var = FlaxEncoder(
            n_input=self.n_input,
            n_latent=self.n_latent,
            n_hidden=self.n_hidden,
            dropout_rate=self.dropout_rate,
        )(x, self.is_training)
        stddev = jnp.sqrt(var) + 1e-4
        z = mean + stddev * jax.random.normal(z_rng, mean.shape)
        dec_input = jnp.concatenate([z, batch], axis=-1)
        rho_unnorm, disp = FlaxDecoder(
            n_input=self.n_input,
            dropout_rate=self.dropout_rate,
            n_hidden=self.n_hidden,
        )(dec_input, self.is_training)
        disp_ = jnp.exp(disp)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho
        nb_logits = jnp.log(mu + 1e-6) - jnp.log(disp_ + 1e-6)
        disp_ = jnp.exp(disp)

        nb = dist.NegativeBinomialLogits(logits=nb_logits, total_count=disp_)

        return VAEOutput(mean, stddev, nb)
