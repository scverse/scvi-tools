import math
from typing import NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS


class Dense(nn.Dense):
    def __init__(self, *args, **kwargs):
        # https://github.com/google/jax/blob/ab15db7d8658825afa9709df46f0dea76688d309/jax/_src/nn/initializers.py#L169
        # sqrt 5 comes from PyTorch
        init = variance_scaling(math.sqrt(5.0), "fan_in", "uniform")
        kwargs.update({"kernel_init": init})
        super().__init__(*args, **kwargs)


class FlaxEncoder(nn.Module):
    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int
    is_training: Optional[bool] = None

    @nn.compact
    def __call__(self, inputs, is_training):

        is_training = nn.merge_param("is_training", self.is_training, is_training)

        inputs_ = jnp.log1p(inputs)

        h = Dense(self.n_hidden)(inputs_)
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        h = nn.Dropout(self.dropout_rate)(h, deterministic=not is_training)
        h = Dense(self.n_hidden)(h)
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        h = nn.Dropout(self.dropout_rate)(h, deterministic=not is_training)

        mean = Dense(self.n_latent)(h)
        log_concentration = Dense(self.n_latent)(h)

        return mean, nn.softplus(log_concentration)


class FlaxDecoder(nn.Module):
    n_input: int
    dropout_rate: float
    n_hidden: int
    is_training: Optional[bool] = None

    @nn.compact
    def __call__(self, z, batch, is_training):
        disp = self.param("disp", lambda rng, shape: jnp.ones(shape), (self.n_input, 1))
        is_training = nn.merge_param("is_training", self.is_training, is_training)
        h = Dense(self.n_hidden)(z)
        h += Dense(self.n_hidden)(batch)

        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        h = nn.Dropout(self.dropout_rate)(h, deterministic=not is_training)
        # skip connection
        h = Dense(self.n_hidden)(jnp.concatenate([h, batch], axis=-1))
        h = nn.LayerNorm(
            use_scale=False,
            use_bias=False,
        )(h)
        h = nn.relu(h)
        h = nn.Dropout(self.dropout_rate)(h, deterministic=not is_training)
        h = Dense(self.n_input)(h)
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
    def __call__(self, array_dict) -> VAEOutput:

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

        qz = dist.Normal(mean, stddev)
        z_rng = self.make_rng("z")
        z = qz.rsample(z_rng)
        rho_unnorm, disp = FlaxDecoder(
            n_input=self.n_input,
            dropout_rate=self.dropout_rate,
            n_hidden=self.n_hidden,
        )(z, batch, self.is_training)
        disp_ = jnp.exp(disp)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho
        nb_logits = jnp.log(mu + 1e-6) - jnp.log(disp_ + 1e-6)
        disp_ = jnp.exp(disp)

        nb = dist.NegativeBinomialLogits(logits=nb_logits, total_count=disp_)

        return VAEOutput(mean, stddev, nb)
