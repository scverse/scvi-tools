import jax
import jax.numpy as jnp

import numpyro.distributions as dist
from scvi.module.base import JAXBaseModuleClass, LossRecorder


class JaxPeakVAE(JAXBaseModuleClass):
    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    n_layers: int = 1
    eps: float = 1e-8

    def setup(self):
        dist.Bernoulli

    def 
