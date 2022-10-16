from functools import partial
from typing import Dict, NamedTuple, Union

import jax
import jax.numpy as jnp

from scvi.module.base import JaxBaseModuleClass, LossRecorder


class _TANGRAM_REGISTRY_KEYS_NT(NamedTuple):
    SC_KEY: str = "X"
    SP_KEY: str = "Y"
    DENSITY_KEY: str = "DENSITY"


TANGRAM_REGISTRY_KEYS = _TANGRAM_REGISTRY_KEYS_NT()

EPS = 1e-8


@partial(jax.jit, static_argnums=(1,))
def _squared_norm(
    x: Union[jnp.ndarray, jax.experimental.sparse.BCOO], axis: int
) -> jnp.ndarray:
    squared = jnp.multiply(x, x)
    res = squared.sum(axis=axis)
    res = res @ jnp.identity(res.shape[0])
    return res.ravel()


@partial(jax.jit, static_argnums=(2,))
def _cosine_similarity(
    x: Union[jnp.ndarray, jax.experimental.sparse.BCOO], y: jnp.ndarray, axis: int
) -> jnp.ndarray:
    """Compute cosine similarity.

    Any weirdness here is because of the sparse matrix.
    """
    x_norm = _squared_norm(x, axis=axis)
    y_norm = _squared_norm(y, axis=axis)
    if axis == 0:
        numerator = jnp.diagonal((y.transpose() @ x).transpose())
    elif axis == 1:
        numerator = jnp.diagonal(x @ y.tranpose())
    else:
        raise ValueError("axis must be 0 or 1")
    denom = jnp.sqrt(x_norm * y_norm)
    return numerator / (denom + EPS)


def _density_criterion(log_y_pred: jnp.ndarray, y_true: jnp.ndarray) -> jnp.ndarray:
    # Kl divergence between the predicted and true distributions
    log_y_true = jnp.log(y_true + EPS)
    return (y_true * (log_y_true - log_y_pred)).sum()


class TangramMapper(JaxBaseModuleClass):
    """Tangram Mapper Model."""

    n_obs_sc: int
    n_obs_sp: int
    lambda_g1: float = 1.0
    lambda_d: float = 0
    lambda_g2: float = 0
    lambda_r: float = 0

    def setup(self):
        """Setup model."""
        self.mapper_unconstrained = self.param(
            "mapper_unconstrained",
            lambda rng, shape: jax.random.normal(rng, shape),
            (self.n_obs_sp, self.n_obs_sc),
        )

    @property
    def required_rngs(self):  # noqa: D102
        return ("params",)

    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray]):
        """Get input for inference."""
        return {}

    def inference(self) -> dict:
        """Run inference model."""
        return {}

    def _get_generative_input(
        self,
        tensors: Dict[str, jnp.ndarray],
        inference_outputs: Dict[str, jnp.ndarray],
    ):
        return {}

    def generative(self) -> dict:
        """No generative model here."""
        return {}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
    ):
        """Compute loss."""
        sp = tensors[TANGRAM_REGISTRY_KEYS.SP_KEY]
        sc = tensors[TANGRAM_REGISTRY_KEYS.SC_KEY]
        mapper = jax.nn.softmax(self.mapper_unconstrained, axis=1)

        if self.lambda_d > 0:
            density = tensors[TANGRAM_REGISTRY_KEYS.DENSITY_KEY]
            d_pred = jnp.log(mapper.sum(axis=0) / mapper.shape[0])
            density_term = self.lambda_d * _density_criterion(d_pred, density)
        else:
            density_term = 0

        g_pred = mapper.transpose() @ sp

        # TODO(adamgayoso): Use these similarities for dense inputs
        # cosine_similarity_0 = jax.vmap(_cosine_similarity_vectors, in_axes=1)
        # cosine_similarity_1 = jax.vmap(_cosine_similarity_vectors, in_axes=0)

        if self.lambda_g1 > 0:
            gv_term = self.lambda_g1 * _cosine_similarity(g_pred, sc, axis=0).mean()
        else:
            gv_term = 0
        if self.lambda_g2 > 0:
            vg_term = self.lambda_g2 * _cosine_similarity(g_pred, sc, axis=1).mean()
        else:
            vg_term = 0

        expression_term = gv_term + vg_term

        regularizer_term = self.lambda_r * (jnp.log(mapper) * mapper).sum()

        total_loss = -expression_term - regularizer_term
        total_loss = total_loss + density_term

        return LossRecorder(total_loss, expression_term, regularizer_term)
