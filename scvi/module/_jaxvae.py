from typing import Dict, NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial


class Dense(nn.Dense):
    def __init__(self, *args, **kwargs):
        # https://github.com/google/jax/blob/ab15db7d8658825afa9709df46f0dea76688d309/jax/_src/nn/initializers.py#L169
        # sqrt 5 comes from PyTorch
        kernel_init = variance_scaling(5.0, "fan_in", "uniform")
        bias_init = variance_scaling(5.0, "fan_in", "uniform", in_axis=-1)
        kwargs.update({"kernel_init": kernel_init, "bias_init": bias_init})
        super().__init__(*args, **kwargs)


class FlaxEncoder(nn.Module):
    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int
    is_training: Optional[bool] = None

    def setup(self):
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_latent)
        self.dense4 = Dense(self.n_latent)

        training = not self.is_training
        self.batchnorm1 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.batchnorm2 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.dropout1 = nn.Dropout(self.dropout_rate, deterministic=training)
        self.dropout2 = nn.Dropout(self.dropout_rate, deterministic=training)

    def __call__(self, x: jnp.ndarray, is_training: bool):

        x_ = jnp.log1p(x)

        h = self.dense1(x_)
        h = self.batchnorm1(h)
        h = nn.relu(h)
        h = self.dropout1(h)
        h = self.dense2(h)
        h = self.batchnorm2(h)
        h = nn.relu(h)
        h = self.dropout2(h)

        mean = self.dense3(h)
        log_var = self.dense4(h)

        return mean, jnp.exp(log_var)


class FlaxDecoder(nn.Module):
    n_input: int
    dropout_rate: float
    n_hidden: int
    is_training: Optional[bool] = None

    def setup(self):
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_hidden)
        self.dense4 = Dense(self.n_input)

        training = not self.is_training
        self.batchnorm1 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.batchnorm2 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.dropout1 = nn.Dropout(self.dropout_rate, deterministic=training)
        self.dropout2 = nn.Dropout(self.dropout_rate, deterministic=training)

    def __call__(self, z: jnp.ndarray, batch: jnp.ndarray, is_training: bool):
        disp = self.param(
            "disp", lambda rng, shape: jax.random.normal(rng, shape), (self.n_input, 1)
        )
        h = self.dense1(z)
        h += self.dense2(batch)

        h = self.batchnorm1(h)
        h = nn.relu(h)
        h = self.dropout1(h)
        # skip connection
        h = self.dense3(jnp.concatenate([h, batch], axis=-1))
        h = self.batchnorm2(h)
        h = nn.relu(h)
        h = self.dropout2(h)
        h = self.dense4(self.n_input)(h)
        return h, disp.ravel()


class VAEOutput(NamedTuple):
    rec_loss: jnp.ndarray
    kl: jnp.ndarray
    px: NegativeBinomial
    qz: dist.Normal


class JaxVAE(nn.Module):
    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    is_training: bool = False
    n_layers: int = 1
    gene_likelihood: str = "nb"
    eps: float = 1e-8

    def setup(self):
        self.encoder = FlaxEncoder(
            n_input=self.n_input,
            n_latent=self.n_latent,
            n_hidden=self.n_hidden,
            dropout_rate=self.dropout_rate,
            is_training=self.is_training,
        )

        self.decoder = FlaxDecoder(
            n_input=self.n_input,
            dropout_rate=0.0,
            n_hidden=self.n_hidden,
            is_training=self.is_training,
        )

    def __call__(self, array_dict: Dict[str, np.ndarray]) -> VAEOutput:

        x = array_dict[REGISTRY_KEYS.X_KEY]
        batch = array_dict[REGISTRY_KEYS.BATCH_KEY]

        # one hot adds an extra dimension
        batch = jax.nn.one_hot(batch, self.n_batch).squeeze(-2)
        mean, var = self.encoder(x)
        stddev = jnp.sqrt(var) + self.eps

        qz = dist.Normal(mean, stddev)
        z_rng = self.make_rng("z")
        z = qz.rsample(z_rng)
        rho_unnorm, disp = self.decoder(z, batch)
        disp_ = jnp.exp(disp)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho

        if self.gene_likelihood == "nb":
            disp_ = jnp.exp(disp)
            px = NegativeBinomial(mean=mu, inverse_dispersion=disp_)
        else:
            px = dist.Poisson(mu)

        rec_loss = -px.log_prob(x).sum(-1)
        kl_div = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)

        return VAEOutput(rec_loss, kl_div, px, qz)
