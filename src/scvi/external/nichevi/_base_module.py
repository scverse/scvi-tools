from __future__ import annotations

from dataclasses import field
from typing import TYPE_CHECKING

import flax

if TYPE_CHECKING:
    from collections.abc import Iterable

    from scvi._types import LossRecord, Tensor


@flax.struct.dataclass
class NicheLossOutput:
    """Loss signature for models.

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
    composition_loss: LossRecord | None = None
    niche_loss: LossRecord | None = None
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
        if self.composition_loss is None:
            object.__setattr__(self, "composition_loss", default)
        if self.niche_loss is None:
            object.__setattr__(self, "niche_loss", default)

        object.__setattr__(self, "reconstruction_loss", self._as_dict("reconstruction_loss"))
        object.__setattr__(self, "kl_local", self._as_dict("kl_local"))
        object.__setattr__(self, "kl_global", self._as_dict("kl_global"))
        object.__setattr__(self, "composition_loss", self._as_dict("composition_loss"))
        object.__setattr__(self, "niche_loss", self._as_dict("niche_loss"))
        object.__setattr__(
            self,
            "reconstruction_loss_sum",
            self.dict_sum(self.reconstruction_loss).sum(),
        )
        object.__setattr__(self, "kl_local_sum", self.dict_sum(self.kl_local).sum())
        object.__setattr__(self, "kl_global_sum", self.dict_sum(self.kl_global))
        object.__setattr__(self, "classification_loss", self.dict_sum(self.classification_loss))
        object.__setattr__(self, "extra_metrics", self.extra_metrics)

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
