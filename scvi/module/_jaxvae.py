from typing import Dict, Optional

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial
from scvi.module.base import JaxBaseModuleClass, LossRecorder


class Dense(nn.Dense):
    def __init__(self, *args, **kwargs):
        # scale set to reimplement pytorch init
        scale = 1 / 3
        kernel_init = variance_scaling(scale, "fan_in", "uniform")
        # bias init can't see input shape so don't include here
        kwargs.update({"kernel_init": kernel_init})
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

    def __call__(self, x: jnp.ndarray):

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
        self.dense4 = Dense(self.n_hidden)
        self.dense5 = Dense(self.n_input)

        training = not self.is_training
        self.batchnorm1 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.batchnorm2 = nn.BatchNorm(
            momentum=0.99, epsilon=0.001, use_running_average=training
        )
        self.dropout1 = nn.Dropout(self.dropout_rate, deterministic=training)
        self.dropout2 = nn.Dropout(self.dropout_rate, deterministic=training)

        self.disp = self.param(
            "disp", lambda rng, shape: jax.random.normal(rng, shape), (self.n_input, 1)
        )

    def __call__(self, z: jnp.ndarray, batch: jnp.ndarray):

        h = self.dense1(z)
        h += self.dense2(batch)

        h = self.batchnorm1(h)
        h = nn.relu(h)
        h = self.dropout1(h)
        h = self.dense3(h)
        # skip connection
        h += self.dense4(batch)
        h = self.batchnorm2(h)
        h = nn.relu(h)
        h = self.dropout2(h)
        h = self.dense5(h)
        return h, self.disp.ravel()


class JaxVAE(JaxBaseModuleClass):
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

    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray]):
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = dict(x=x)
        return input_dict

    def inference(self, x: jnp.ndarray, n_samples: int = 1) -> dict:
        mean, var = self.encoder(x)
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

        # one hot adds an extra dimension
        batch = jax.nn.one_hot(batch_index, self.n_batch).squeeze(-2)
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

        return dict(px=px, rho=rho)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        x = tensors[REGISTRY_KEYS.X_KEY]
        px = generative_outputs["px"]
        qz = inference_outputs["qz"]
        reconst_loss = -px.log_prob(x).sum(-1)
        kl_divergence_z = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)

        kl_local_for_warmup = kl_divergence_z
        weighted_kl_local = kl_weight * kl_local_for_warmup

        loss = jnp.mean(reconst_loss + weighted_kl_local)

        kl_local = kl_divergence_z
        return LossRecorder(loss, reconst_loss, kl_local)
