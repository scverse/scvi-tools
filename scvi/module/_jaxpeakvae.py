from typing import Dict

import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial
from scvi.module.base import JaxBaseModuleClass, LossRecorder

from ._jaxvae import FlaxDecoder, FlaxEncoder


class JaxPEAKVAE(JaxBaseModuleClass):
    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    n_layers: int = 1
    eps: float = 1e-8

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

    @property
    def required_rngs(self):
        return ("params", "dropout", "z")
    
    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray]):
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = dict(x=x)
        return input_dict
    
    def inference(self, x: jnp.ndarray, n_samples: int = 1) -> dict:
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
        rho_unnorm, disp = self.decoder(z, batch, training=self.training)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho

        px = dist.Bernoulli(mu)

        return dict(px=px)
    
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

        reconstruction_loss = -px.log_prob(x).sum(-1)

        kl_div = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)
        loss = jnp.mean(reconstruction_loss + kl_weight * kl_div)
        

        return LossRecorder(loss, reconstruction_loss, kl_div)

    
