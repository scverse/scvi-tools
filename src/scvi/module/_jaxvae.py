from collections.abc import Iterable
from dataclasses import field

import flax
import jax
import jax.numpy as jnp
import numpyro.distributions as dist
from flax import linen as nn
from flax.linen.initializers import variance_scaling

from scvi import REGISTRY_KEYS
from scvi._types import LossRecord, Tensor
from scvi.distributions import JaxNegativeBinomialMeanDisp as NegativeBinomial
from scvi.module.base import JaxBaseModuleClass, flax_configure


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


class Dense(nn.Dense):
    """Jax dense layer."""

    def __init__(self, *args, **kwargs):
        # scale set to reimplement pytorch init
        scale = 1 / 3
        kernel_init = variance_scaling(scale, "fan_in", "uniform")
        # bias init can't see input shape so don't include here
        kwargs.update({"kernel_init": kernel_init})
        super().__init__(*args, **kwargs)


class FlaxEncoder(nn.Module):
    """Encoder for Jax VAE."""

    n_input: int
    n_latent: int
    n_hidden: int
    dropout_rate: int
    training: bool | None = None

    def setup(self):
        """Setup encoder."""
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_latent)
        self.dense4 = Dense(self.n_latent)

        self.batchnorm1 = nn.BatchNorm(momentum=0.9)
        self.batchnorm2 = nn.BatchNorm(momentum=0.9)
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

    def __call__(self, x: jnp.ndarray, training: bool | None = None):
        """Forward pass."""
        training = nn.merge_param("training", self.training, training)
        is_eval = not training

        x_ = jnp.log1p(x)

        h = self.dense1(x_)
        h = self.batchnorm1(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense2(h)
        h = self.batchnorm2(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout2(h, deterministic=is_eval)

        mean = self.dense3(h)
        log_var = self.dense4(h)

        return mean, jnp.exp(log_var)


class FlaxDecoder(nn.Module):
    """Decoder for Jax VAE."""

    n_input: int
    dropout_rate: float
    n_hidden: int
    training: bool | None = None

    def setup(self):
        """Setup decoder."""
        self.dense1 = Dense(self.n_hidden)
        self.dense2 = Dense(self.n_hidden)
        self.dense3 = Dense(self.n_hidden)
        self.dense4 = Dense(self.n_hidden)
        self.dense5 = Dense(self.n_input)

        self.batchnorm1 = nn.BatchNorm(momentum=0.9)
        self.batchnorm2 = nn.BatchNorm(momentum=0.9)
        self.dropout1 = nn.Dropout(self.dropout_rate)
        self.dropout2 = nn.Dropout(self.dropout_rate)

        self.disp = self.param(
            "disp", lambda rng, shape: jax.random.normal(rng, shape), (self.n_input, 1)
        )

    def __call__(self, z: jnp.ndarray, batch: jnp.ndarray, training: bool | None = None):
        """Forward pass."""
        # TODO(adamgayoso): Test this
        training = nn.merge_param("training", self.training, training)
        is_eval = not training

        h = self.dense1(z)
        h += self.dense2(batch)

        h = self.batchnorm1(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout1(h, deterministic=is_eval)
        h = self.dense3(h)
        # skip connection
        h += self.dense4(batch)
        h = self.batchnorm2(h, use_running_average=is_eval)
        h = nn.relu(h)
        h = self.dropout2(h, deterministic=is_eval)
        h = self.dense5(h)
        return h, self.disp.ravel()


@flax_configure
class JaxVAE(JaxBaseModuleClass):
    """Variational autoencoder model."""

    n_input: int
    n_batch: int
    n_hidden: int = 128
    n_latent: int = 30
    dropout_rate: float = 0.0
    n_layers: int = 1
    gene_likelihood: str = "nb"
    eps: float = 1e-8
    training: bool = True

    def setup(self):
        """Setup model."""
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

    def _get_inference_input(self, tensors: dict[str, jnp.ndarray]):
        """Get input for inference."""
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = {"x": x}
        return input_dict

    def inference(self, x: jnp.ndarray, n_samples: int = 1) -> dict:
        """Run inference model."""
        mean, var = self.encoder(x, training=self.training)
        stddev = jnp.sqrt(var) + self.eps

        qz = dist.Normal(mean, stddev)
        z_rng = self.make_rng("z")
        sample_shape = () if n_samples == 1 else (n_samples,)
        z = qz.rsample(z_rng, sample_shape=sample_shape)

        return {"qz": qz, "z": z}

    def _get_generative_input(
        self,
        tensors: dict[str, jnp.ndarray],
        inference_outputs: dict[str, jnp.ndarray],
    ):
        """Get input for generative model."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        z = inference_outputs["z"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            "x": x,
            "z": z,
            "batch_index": batch_index,
        }
        return input_dict

    def generative(self, x, z, batch_index) -> dict:
        """Run generative model."""
        # one hot adds an extra dimension
        batch = jax.nn.one_hot(batch_index, self.n_batch).squeeze(-2)
        rho_unnorm, disp = self.decoder(z, batch, training=self.training)
        disp_ = jnp.exp(disp)
        rho = jax.nn.softmax(rho_unnorm, axis=-1)
        total_count = x.sum(-1)[:, jnp.newaxis]
        mu = total_count * rho

        if self.gene_likelihood == "nb":
            disp_ = jnp.exp(disp)
            px = NegativeBinomial(mean=mu, inverse_dispersion=disp_)
        else:
            px = dist.Poisson(mu)

        return {"px": px, "rho": rho}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Compute loss."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        px = generative_outputs["px"]
        qz = inference_outputs["qz"]
        reconst_loss = -px.log_prob(x).sum(-1)
        kl_divergence_z = dist.kl_divergence(qz, dist.Normal(0, 1)).sum(-1)

        kl_local_for_warmup = kl_divergence_z
        weighted_kl_local = kl_weight * kl_local_for_warmup

        loss = jnp.mean(reconst_loss + weighted_kl_local)

        kl_local = kl_divergence_z
        return LossOutput(loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_local)
