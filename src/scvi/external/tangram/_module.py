from collections.abc import Iterable
from dataclasses import field
from typing import NamedTuple

import flax
import jax
import jax.numpy as jnp

from scvi._types import LossRecord, Tensor
from scvi.module.base._base_module import JaxBaseModuleClass
from scvi.module.base._decorators import flax_configure


@flax.struct.dataclass
class LossOutput:
    """Loss signature for Jax models.

    This class provides an organized way to record the model loss, as well as
    the components of the ELBO. This may also be used in MLE, MAP, EM methods.
    The loss is used for backpropagation during inference. The other parameters
    are used for logging/early stopping during inference.

    Parameters
    ----------
    loss
        Tensor with loss for minibatch. Should be one dimensional with one value.
        Note that loss should be in an array/tensor and not a float.
    reconstruction_loss
        Reconstruction loss for each observation in the minibatch. If a tensor, converted to
        a dictionary with key "reconstruction_loss" and value as tensor.
    kl_local
        KL divergence associated with each observation in the minibatch. If a tensor, converted to
        a dictionary with key "kl_local" and value as tensor.
    kl_global
        Global KL divergence term. Should be one dimensional with one value. If a tensor, converted
        to a dictionary with key "kl_global" and value as tensor.
    classification_loss
        Classification loss.
    logits
        Logits for classification.
    true_labels
        True labels for classification.
    extra_metrics
        Additional metrics can be passed as arrays/tensors or dictionaries of
        arrays/tensors.
    n_obs_minibatch
        Number of observations in the minibatch. If None, will be inferred from
        the shape of the reconstruction_loss tensor.


    Examples
    --------
    >>> loss_output = LossOutput(
    ...     loss=loss,
    ...     reconstruction_loss=reconstruction_loss,
    ...     kl_local=kl_local,
    ...     extra_metrics={"x": scalar_tensor_x, "y": scalar_tensor_y},
    ... )
    """

    loss: LossRecord
    reconstruction_loss: LossRecord | None = None
    kl_local: LossRecord | None = None
    kl_global: LossRecord | None = None
    classification_loss: LossRecord | None = None
    logits: Tensor | None = None
    true_labels: Tensor | None = None
    extra_metrics: dict[str, Tensor] | None = field(default_factory=dict)
    n_obs_minibatch: int | None = None
    reconstruction_loss_sum: Tensor = field(default=None)
    kl_local_sum: Tensor = field(default=None)
    kl_global_sum: Tensor = field(default=None)

    def __post_init__(self):
        object.__setattr__(self, "loss", self.dict_sum(self.loss))

        if self.n_obs_minibatch is None and self.reconstruction_loss is None:
            raise ValueError("Must provide either n_obs_minibatch or reconstruction_loss")

        default = 0 * self.loss
        if self.reconstruction_loss is None:
            object.__setattr__(self, "reconstruction_loss", default)
        if self.kl_local is None:
            object.__setattr__(self, "kl_local", default)
        if self.kl_global is None:
            object.__setattr__(self, "kl_global", default)

        object.__setattr__(self, "reconstruction_loss", self._as_dict("reconstruction_loss"))
        object.__setattr__(self, "kl_local", self._as_dict("kl_local"))
        object.__setattr__(self, "kl_global", self._as_dict("kl_global"))
        object.__setattr__(
            self,
            "reconstruction_loss_sum",
            self.dict_sum(self.reconstruction_loss).sum(),
        )
        object.__setattr__(self, "kl_local_sum", self.dict_sum(self.kl_local).sum())
        object.__setattr__(self, "kl_global_sum", self.dict_sum(self.kl_global))

        if self.reconstruction_loss is not None and self.n_obs_minibatch is None:
            rec_loss = self.reconstruction_loss
            object.__setattr__(self, "n_obs_minibatch", list(rec_loss.values())[0].shape[0])

        if self.classification_loss is not None and (
            self.logits is None or self.true_labels is None
        ):
            raise ValueError(
                "Must provide `logits` and `true_labels` if `classification_loss` is provided."
            )

    @staticmethod
    def dict_sum(dictionary: dict[str, Tensor] | Tensor):
        """Sum over elements of a dictionary."""
        if isinstance(dictionary, dict):
            return sum(dictionary.values())
        else:
            return dictionary

    @property
    def extra_metrics_keys(self) -> Iterable[str]:
        """Keys for extra metrics."""
        return self.extra_metrics.keys()

    def _as_dict(self, attr_name: str):
        attr = getattr(self, attr_name)
        if isinstance(attr, dict):
            return attr
        else:
            return {attr_name: attr}


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
    target_count: int | None = None
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
    def required_rngs(self):
        return ("params",)

    def _get_inference_input(self, tensors: dict[str, jnp.ndarray]):
        """Get input for inference."""
        return {}

    def inference(self) -> dict:
        """Run inference model."""
        return {}

    def _get_generative_input(
        self,
        tensors: dict[str, jnp.ndarray],
        inference_outputs: dict[str, jnp.ndarray],
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
