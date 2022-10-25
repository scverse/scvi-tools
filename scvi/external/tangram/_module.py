from typing import Dict, NamedTuple, Optional

import chex
import jax
import jax.numpy as jnp

from scvi.module.base import JaxBaseModuleClass, LossOutput, flax_configure


class _TANGRAM_REGISTRY_KEYS_NT(NamedTuple):
    SC_KEY: str = "X"
    SP_KEY: str = "Y"
    DENSITY_KEY: str = "DENSITY"


TANGRAM_REGISTRY_KEYS = _TANGRAM_REGISTRY_KEYS_NT()

EPS = 1e-8


def _cosine_similarity_vectors(x: jnp.ndarray, y: jnp.ndarray) -> jnp.ndarray:
    return jnp.dot(x, y) / (jnp.maximum(jnp.linalg.norm(x) * jnp.linalg.norm(y), EPS))


def _density_criterion(log_y_pred: jnp.ndarray, y_true: jnp.ndarray) -> jnp.ndarray:
    # Kl divergence between the predicted and true distributions
    log_y_true = jnp.log(y_true + EPS)
    return (y_true * (log_y_true - log_y_pred)).sum()


@flax_configure
class TangramMapper(JaxBaseModuleClass):
    """Tangram Mapper Model."""

    n_obs_sc: int
    n_obs_sp: int
    lambda_g1: float = 1.0
    lambda_d: float = 0.0
    lambda_g2: float = 0.0
    lambda_r: float = 0.0
    lambda_count: float = 1.0
    lambda_f_reg: float = 1.0
    constrained: bool = False
    target_count: Optional[int] = None
    training: bool = True

    def setup(self):
        """Setup model."""
        self.mapper_unconstrained = self.param(
            "mapper_unconstrained",
            lambda rng, shape: jax.random.normal(rng, shape),
            (self.n_obs_sc, self.n_obs_sp),
        )

        if self.constrained:
            self.filter_unconstrained = self.param(
                "filter_unconstrained",
                lambda rng, shape: jax.random.normal(rng, shape),
                (self.n_obs_sc, 1),
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

        if self.constrained:
            filter = jax.nn.sigmoid(self.filter_unconstrained)
            mapper_filtered = mapper * filter

        if self.lambda_d > 0:
            density = tensors[TANGRAM_REGISTRY_KEYS.DENSITY_KEY].ravel()
            if self.constrained:
                d_pred = jnp.log(mapper_filtered.sum(axis=0) / (filter.sum()))
            else:
                d_pred = jnp.log(mapper.sum(axis=0) / mapper.shape[0])
            density_term = self.lambda_d * _density_criterion(d_pred, density)
        else:
            density_term = 0

        if self.constrained:
            sc = sc * filter

        g_pred = mapper.transpose() @ sc
        chex.assert_equal_shape([sp, g_pred])

        # Expression term
        if self.lambda_g1 > 0:
            cosine_similarity_0 = jax.vmap(_cosine_similarity_vectors, in_axes=1)
            gv_term = self.lambda_g1 * cosine_similarity_0(sp, g_pred).mean()
        else:
            gv_term = 0
        if self.lambda_g2 > 0:
            cosine_similarity_1 = jax.vmap(_cosine_similarity_vectors, in_axes=0)
            vg_term = self.lambda_g1 * cosine_similarity_1(sp, g_pred).mean()
            vg_term = self.lambda_g2 * vg_term
        else:
            vg_term = 0

        expression_term = gv_term + vg_term

        # Regularization terms
        if self.lambda_r > 0:
            regularizer_term = self.lambda_r * (jnp.log(mapper) * mapper).sum()
        else:
            regularizer_term = 0

        if self.lambda_count > 0 and self.constrained:
            if self.target_count is None:
                raise ValueError("target_count must be set if in constrained mode.")
            count_term = self.lambda_count * jnp.abs(filter.sum() - self.target_count)
        else:
            count_term = 0

        if self.lambda_f_reg > 0 and self.constrained:
            f_reg_t = filter - jnp.square(filter)
            f_reg = self.lambda_f_reg * f_reg_t.sum()
        else:
            f_reg = 0

        # Total loss
        total_loss = -expression_term - regularizer_term + count_term + f_reg
        total_loss = total_loss + density_term

        return LossOutput(
            loss=total_loss,
            n_obs_minibatch=sp.shape[0],
            extra_metrics={
                "expression_term": expression_term,
                "regularizer_term": regularizer_term,
            },
        )
