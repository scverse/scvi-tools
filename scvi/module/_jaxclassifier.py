from typing import Optional

import jax
import jax.numpy as jnp
from flax import linen as nn
from scvi.module._jaxvae import Dense


class FlaxClassifier(nn.Module):
    """Basic fully-connected NN classifier Flax module.

    Parameters
    ----------
    n_labels
        Numput of outputs dimensions
    n_hidden
        Number of hidden nodes in hidden layer
    n_layers
        Number of hidden layers
    logits
        Return logits or not
    """

    n_labels: int
    n_hidden: int = 32
    n_layers: int = 2
    logits: bool = True
    training: Optional[bool] = None

    @nn.compact
    def __call__(self, x: jnp.ndarray, training: Optional[bool] = None):
        training = nn.merge_param("training", self.training, training)

        h = x
        for _ in range(self.n_layers - 1):
            h = Dense(self.n_hidden)(h)
            h = nn.relu(h)
        out = Dense(self.n_labels)(h)

        if not self.logits:
            out = jax.nn.softmax(out, axis=-1)

        return out
