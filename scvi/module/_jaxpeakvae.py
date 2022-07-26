from typing import Dict

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
import optax
from flax import linen as nn

from scvi import REGISTRY_KEYS
from scvi.module.base import JaxBaseModuleClass, LossRecorder

from ._jaxvae import Dense


class FlaxEncoder(nn.Module):
    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int

    def setup(self):
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_latent)
        self.dense4 = Dense(self.n_latent)

        self.layernorm1 = nn.LayerNorm()
        self.layernorm2 = nn.LayerNorm()
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

    def __call__(self, x: jnp.ndarray, training: bool = False):
        is_eval = not training

        x_ = jnp.log1p(x)

        h = self.dense1(x_)
        h = self.layernorm1(h)
        h = nn.leaky_relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense2(h)
        h = self.layernorm2(h)
        h = nn.leaky_relu(h)
        h = self.dropout2(h, deterministic=is_eval)

        mean = self.dense3(h)
        log_var = self.dense4(h)

        return mean, jnp.exp(log_var)


class FlaxDecoder(nn.Module):
    n_input: int
    dropout_rate: float
    n_hidden: int

    def setup(self):
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_hidden)
        self.dense4 = Dense(self.n_hidden)
        self.dense5 = Dense(self.n_input)

        self.batchnorm1 = nn.BatchNorm()
        self.batchnorm2 = nn.BatchNorm()
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

    def __call__(self, z: jnp.ndarray, batch: jnp.ndarray, training: bool = False):
        is_eval = not training

        h = self.dense1(z)
        h += self.dense2(batch)

        h = self.batchnorm1(h, use_running_average=is_eval)
        h = nn.leaky_relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense3(h)
        # skip connection
        h += self.dense4(batch)
        h = self.batchnorm2(h, use_running_average=is_eval)
        h = nn.leaky_relu(h)
        h = self.dropout2(h, deterministic=is_eval)
        h = self.dense5(h)
        h = nn.sigmoid(h)

        return h


class JaxPEAKVAE(JaxBaseModuleClass):
    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    n_layers: int = 1
    eps: float = 1e-8
    region_factors: bool = True

    def setup(self):
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

        self.d_encoder = FlaxDecoder(
            n_input=self.n_latent,
            dropout_rate=0.0,
            n_hidden=self.n_hidden,
        )

        if self.region_factors:
            self.rf = self.param("rf", nn.initializers.zeros, self.n_input)
        else:
            self.rf = 1

    @property
    def required_rngs(self):
        return ("params", "dropout", "z")

    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray]):
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = dict(x=x)
        return input_dict

    def inference(self, x: jnp.ndarray, batch_index, n_samples: int = 1) -> dict:
        mean, var = self.encoder(x, training=self.training)
        stddev = jnp.sqrt(var) + self.eps

        batch = jax.nn.one_hot(batch_index, self.n_batch).squeeze(-2)
        d = self.d_encoder(x, batch, training=self.training)

        qz = dist.Normal(mean, stddev)
        z_rng = self.make_rng("z")
        sample_shape = () if n_samples == 1 else (n_samples,)
        z = qz.rsample(z_rng, sample_shape=sample_shape)

        return dict(qz=qz, z=z, d=d)

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
        batch = jax.nn.one_hot(batch_index, self.n_batch).squeeze(-2)
        rho = self.decoder(z, batch, training=self.training)
        # total_count = x.sum(-1)[:, jnp.newaxis]
        # mu = total_count * rho

        # px = dist.Bernoulli(mu) # numpyro distribution object

        return dict(px=rho)
    
    def bceloss(logits, labels, eps=1e-8):
        return -labels * jnp.log(logits + eps) - (1.0 - labels) * jnp.log(1 - logits + eps)
    
    def get_reconstruction_loss(self, p, d, f, x):
        return self.bceloss(p * d * f, (x > 0).float()).sum(dim=-1)

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
        d = inference_outputs["d"]

        print('777777777777777')
        print(self.rf * px)

        reconstruction_loss = -px.log_prob(x).sum(-1)

        kl_div = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)
        loss = jnp.mean(reconstruction_loss + kl_weight * kl_div)

        

        return LossRecorder(loss, reconstruction_loss, kl_div)
