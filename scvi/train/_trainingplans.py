from collections import OrderedDict
from functools import partial
from inspect import signature
from typing import Callable, Dict, Iterable, Literal, Optional, Union

import jax
import jax.numpy as jnp
import numpy as np
import optax
import pyro
import pytorch_lightning as pl
import torch
from pyro.nn import PyroModule
from torch.optim.lr_scheduler import ReduceLROnPlateau

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable, TunableMixin
from scvi.module import Classifier
from scvi.module.base import (
    BaseModuleClass,
    JaxBaseModuleClass,
    LossOutput,
    LossRecorder,
    PyroBaseModuleClass,
    TrainStateWithState,
)
from scvi.nn import one_hot

from ._metrics import ElboMetric

JaxOptimizerCreator = Callable[[], optax.GradientTransformation]
TorchOptimizerCreator = Callable[[Iterable[torch.Tensor]], torch.optim.Optimizer]


def _compute_kl_weight(
    epoch: int,
    step: int,
    n_epochs_kl_warmup: Optional[int],
    n_steps_kl_warmup: Optional[int],
    max_kl_weight: float = 1.0,
    min_kl_weight: float = 0.0,
) -> float:
    """
    Computes the kl weight for the current step or epoch.

    If both `n_epochs_kl_warmup` and `n_steps_kl_warmup` are None `max_kl_weight` is returned.

    Parameters
    ----------
    epoch
        Current epoch.
    step
        Current step.
    n_epochs_kl_warmup
        Number of training epochs to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`
    max_kl_weight
        Maximum scaling factor on KL divergence during training.
    min_kl_weight
        Minimum scaling factor on KL divergence during training.
    """
    if min_kl_weight > max_kl_weight:
        raise ValueError(
            f"min_kl_weight={min_kl_weight} is larger than max_kl_weight={max_kl_weight}."
        )

    slope = max_kl_weight - min_kl_weight
    if n_epochs_kl_warmup:
        if epoch < n_epochs_kl_warmup:
            return slope * (epoch / n_epochs_kl_warmup) + min_kl_weight
    elif n_steps_kl_warmup:
        if step < n_steps_kl_warmup:
            return slope * (step / n_steps_kl_warmup) + min_kl_weight
    return max_kl_weight


class TrainingPlan(TunableMixin, pl.LightningModule):
    """
    Lightning module task to train scvi-tools modules.

    The training plan is a PyTorch Lightning Module that is initialized
    with a scvi-tools module object. It configures the optimizers, defines
    the training step and validation step, and computes metrics to be recorded
    during training. The training step and validation step are functions that
    take data, run it through the model and return the loss, which will then
    be used to optimize the model parameters in the Trainer. Overall, custom
    training plans can be used to develop complex inference schemes on top of
    modules.

    The following developer tutorial will familiarize you more with training plans
    and how to use them: :doc:`/tutorials/notebooks/model_user_guide`.

    Parameters
    ----------
    module
        A module instance from class ``BaseModuleClass``.
    optimizer
        One of "Adam" (:class:`~torch.optim.Adam`), "AdamW" (:class:`~torch.optim.AdamW`),
        or "Custom", which requires a custom optimizer creator callable to be passed via
        `optimizer_creator`.
    optimizer_creator
        A callable taking in parameters and returning a :class:`~torch.optim.Optimizer`.
        This allows using any PyTorch optimizer with custom hyperparameters.
    lr
        Learning rate used for optimization, when `optimizer_creator` is None.
    weight_decay
        Weight decay used in optimization, when `optimizer_creator` is None.
    eps
        eps used for optimization, when `optimizer_creator` is None.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`. Only activated when `n_epochs_kl_warmup` is
        set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from `min_kl_weight` to
        `max_kl_weight`. Overrides `n_steps_kl_warmup` when both are not `None`.
    reduce_lr_on_plateau
        Whether to monitor validation loss and reduce learning rate when validation set
        `lr_scheduler_metric` plateaus.
    lr_factor
        Factor to reduce learning rate.
    lr_patience
        Number of epochs with no improvement after which learning rate will be reduced.
    lr_threshold
        Threshold for measuring the new optimum.
    lr_scheduler_metric
        Which metric to track for learning rate reduction.
    lr_min
        Minimum learning rate allowed.
    max_kl_weight
        Maximum scaling factor on KL divergence during training.
    min_kl_weight
        Minimum scaling factor on KL divergence during training.
    **loss_kwargs
        Keyword args to pass to the loss method of the `module`.
        `kl_weight` should not be passed here and is handled automatically.
    """

    def __init__(
        self,
        module: BaseModuleClass,
        *,
        optimizer: Tunable[Literal["Adam", "AdamW", "Custom"]] = "Adam",
        optimizer_creator: Optional[TorchOptimizerCreator] = None,
        lr: Tunable[float] = 1e-3,
        weight_decay: Tunable[float] = 1e-6,
        eps: Tunable[float] = 0.01,
        n_steps_kl_warmup: Tunable[int] = None,
        n_epochs_kl_warmup: Tunable[int] = 400,
        reduce_lr_on_plateau: Tunable[bool] = False,
        lr_factor: Tunable[float] = 0.6,
        lr_patience: Tunable[int] = 30,
        lr_threshold: Tunable[float] = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        lr_min: Tunable[float] = 0,
        max_kl_weight: Tunable[float] = 1.0,
        min_kl_weight: Tunable[float] = 0.0,
        **loss_kwargs,
    ):
        super().__init__()
        self.module = module
        self.lr = lr
        self.weight_decay = weight_decay
        self.eps = eps
        self.optimizer_name = optimizer
        self.n_steps_kl_warmup = n_steps_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.reduce_lr_on_plateau = reduce_lr_on_plateau
        self.lr_factor = lr_factor
        self.lr_patience = lr_patience
        self.lr_scheduler_metric = lr_scheduler_metric
        self.lr_threshold = lr_threshold
        self.lr_min = lr_min
        self.loss_kwargs = loss_kwargs
        self.min_kl_weight = min_kl_weight
        self.max_kl_weight = max_kl_weight
        self.optimizer_creator = optimizer_creator

        if self.optimizer_name == "Custom" and self.optimizer_creator is None:
            raise ValueError(
                "If optimizer is 'Custom', `optimizer_creator` must be provided."
            )

        self._n_obs_training = None
        self._n_obs_validation = None

        # automatic handling of kl weight
        self._loss_args = set(signature(self.module.loss).parameters.keys())
        if "kl_weight" in self._loss_args:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})

        self.initialize_train_metrics()
        self.initialize_val_metrics()

    @staticmethod
    def _create_elbo_metric_components(mode: str, n_total: Optional[int] = None):
        """Initialize ELBO metric and the metric collection."""
        rec_loss = ElboMetric("reconstruction_loss", mode, "obs")
        kl_local = ElboMetric("kl_local", mode, "obs")
        kl_global = ElboMetric("kl_global", mode, "obs")
        # n_total can be 0 if there is no validation set, this won't ever be used
        # in that case anyway
        n = 1 if n_total is None or n_total < 1 else n_total
        elbo = rec_loss + kl_local + (1 / n) * kl_global
        elbo.name = f"elbo_{mode}"
        collection = OrderedDict(
            [(metric.name, metric) for metric in [elbo, rec_loss, kl_local, kl_global]]
        )
        return elbo, rec_loss, kl_local, kl_global, collection

    def initialize_train_metrics(self):
        """Initialize train related metrics."""
        (
            self.elbo_train,
            self.rec_loss_train,
            self.kl_local_train,
            self.kl_global_train,
            self.train_metrics,
        ) = self._create_elbo_metric_components(
            mode="train", n_total=self.n_obs_training
        )
        self.elbo_train.reset()

    def initialize_val_metrics(self):
        """Initialize val related metrics."""
        (
            self.elbo_val,
            self.rec_loss_val,
            self.kl_local_val,
            self.kl_global_val,
            self.val_metrics,
        ) = self._create_elbo_metric_components(
            mode="validation", n_total=self.n_obs_validation
        )
        self.elbo_val.reset()

    @property
    def n_obs_training(self):
        """
        Number of observations in the training set.

        This will update the loss kwargs for loss rescaling.

        Notes
        -----
        This can get set after initialization
        """
        return self._n_obs_training

    @n_obs_training.setter
    def n_obs_training(self, n_obs: int):
        if "n_obs" in self._loss_args:
            self.loss_kwargs.update({"n_obs": n_obs})
        self._n_obs_training = n_obs
        self.initialize_train_metrics()

    @property
    def n_obs_validation(self):
        """
        Number of observations in the validation set.

        This will update the loss kwargs for loss rescaling.

        Notes
        -----
        This can get set after initialization
        """
        return self._n_obs_validation

    @n_obs_validation.setter
    def n_obs_validation(self, n_obs: int):
        self._n_obs_validation = n_obs
        self.initialize_val_metrics()

    def forward(self, *args, **kwargs):
        """Passthrough to the module's forward method."""
        return self.module(*args, **kwargs)

    @torch.inference_mode()
    def compute_and_log_metrics(
        self,
        loss_recorder: Union[LossRecorder, LossOutput],
        metrics: Dict[str, ElboMetric],
        mode: str,
    ):
        """
        Computes and logs metrics.

        Parameters
        ----------
        loss_recorder
            LossRecorder object from scvi-tools module
        metric_attr_name
            The name of the torch metric object to use
        mode
            Postfix string to add to the metric name of
            extra metrics
        """
        if isinstance(loss_recorder, LossRecorder):
            loss_output = loss_recorder._loss_output
        else:
            loss_output = loss_recorder
        rec_loss = loss_output.reconstruction_loss_sum
        n_obs_minibatch = loss_output.n_obs_minibatch
        kl_local = loss_output.kl_local_sum
        kl_global = loss_output.kl_global_sum

        # Use the torchmetric object for the ELBO
        # We only need to update the ELBO metric
        # As it's defined as a sum of the other metrics
        metrics[f"elbo_{mode}"].update(
            reconstruction_loss=rec_loss,
            kl_local=kl_local,
            kl_global=kl_global,
            n_obs_minibatch=n_obs_minibatch,
        )
        # pytorch lightning handles everything with the torchmetric object
        self.log_dict(
            metrics,
            on_step=False,
            on_epoch=True,
            batch_size=n_obs_minibatch,
        )

        # accumlate extra metrics passed to loss recorder
        for key in loss_output.extra_metrics_keys:
            met = loss_output.extra_metrics[key]
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{key}_{mode}",
                met,
                on_step=False,
                on_epoch=True,
                batch_size=n_obs_minibatch,
            )

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        """Training step for the model."""
        if "kl_weight" in self.loss_kwargs:
            kl_weight = self.kl_weight
            self.loss_kwargs.update({"kl_weight": kl_weight})
            self.log("kl_weight", kl_weight, on_step=True, on_epoch=False)
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("train_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")
        return scvi_loss.loss

    def validation_step(self, batch, batch_idx):
        """Validation step for the model."""
        # loss kwargs here contains `n_obs` equal to n_training_obs
        # so when relevant, the actual loss value is rescaled to number
        # of training examples
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("validation_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.val_metrics, "validation")

    def _optimizer_creator_fn(
        self, optimizer_cls: Union[torch.optim.Adam, torch.optim.AdamW]
    ):
        """
        Create optimizer for the model.

        This type of function can be passed as the `optimizer_creator`
        """
        return lambda params: optimizer_cls(
            params, lr=self.lr, eps=self.eps, weight_decay=self.weight_decay
        )

    def get_optimizer_creator(self):
        """Get optimizer creator for the model."""
        if self.optimizer_name == "Adam":
            optim_creator = self._optimizer_creator_fn(torch.optim.Adam)
        elif self.optimizer_name == "AdamW":
            optim_creator = self._optimizer_creator_fn(torch.optim.AdamW)
        elif self.optimizer_name == "Custom":
            optim_creator = self.optimizer_creator
        else:
            raise ValueError("Optimizer not understood.")

        return optim_creator

    def configure_optimizers(self):
        """Configure optimizers for the model."""
        params = filter(lambda p: p.requires_grad, self.module.parameters())
        optimizer = self.get_optimizer_creator()(params)
        config = {"optimizer": optimizer}
        if self.reduce_lr_on_plateau:
            scheduler = ReduceLROnPlateau(
                optimizer,
                patience=self.lr_patience,
                factor=self.lr_factor,
                threshold=self.lr_threshold,
                min_lr=self.lr_min,
                threshold_mode="abs",
                verbose=True,
            )
            config.update(
                {
                    "lr_scheduler": scheduler,
                    "monitor": self.lr_scheduler_metric,
                },
            )
        return config

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        return _compute_kl_weight(
            self.current_epoch,
            self.global_step,
            self.n_epochs_kl_warmup,
            self.n_steps_kl_warmup,
            self.max_kl_weight,
            self.min_kl_weight,
        )


class AdversarialTrainingPlan(TrainingPlan):
    """
    Train vaes with adversarial loss option to encourage latent space mixing.

    Parameters
    ----------
    module
        A module instance from class ``BaseModuleClass``.
    optimizer
        One of "Adam" (:class:`~torch.optim.Adam`), "AdamW" (:class:`~torch.optim.AdamW`),
        or "Custom", which requires a custom optimizer creator callable to be passed via
        `optimizer_creator`.
    optimizer_creator
        A callable taking in parameters and returning a :class:`~torch.optim.Optimizer`.
        This allows using any PyTorch optimizer with custom hyperparameters.
    lr
        Learning rate used for optimization, when `optimizer_creator` is None.
    weight_decay
        Weight decay used in optimization, when `optimizer_creator` is None.
    eps
        eps used for optimization, when `optimizer_creator` is None.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
        Only activated when `n_epochs_kl_warmup` is set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from 0 to 1.
        Overrides `n_steps_kl_warmup` when both are not `None`.
    reduce_lr_on_plateau
        Whether to monitor validation loss and reduce learning rate when validation set
        `lr_scheduler_metric` plateaus.
    lr_factor
        Factor to reduce learning rate.
    lr_patience
        Number of epochs with no improvement after which learning rate will be reduced.
    lr_threshold
        Threshold for measuring the new optimum.
    lr_scheduler_metric
        Which metric to track for learning rate reduction.
    lr_min
        Minimum learning rate allowed
    adversarial_classifier
        Whether to use adversarial classifier in the latent space
    scale_adversarial_loss
        Scaling factor on the adversarial components of the loss.
        By default, adversarial loss is scaled from 1 to 0 following opposite of
        kl warmup.
    **loss_kwargs
        Keyword args to pass to the loss method of the `module`.
        `kl_weight` should not be passed here and is handled automatically.
    """

    def __init__(
        self,
        module: BaseModuleClass,
        *,
        optimizer: Tunable[Literal["Adam", "AdamW", "Custom"]] = "Adam",
        optimizer_creator: Optional[TorchOptimizerCreator] = None,
        lr: Tunable[float] = 1e-3,
        weight_decay: Tunable[float] = 1e-6,
        n_steps_kl_warmup: Tunable[int] = None,
        n_epochs_kl_warmup: Tunable[int] = 400,
        reduce_lr_on_plateau: Tunable[bool] = False,
        lr_factor: Tunable[float] = 0.6,
        lr_patience: Tunable[int] = 30,
        lr_threshold: Tunable[float] = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        lr_min: float = 0,
        adversarial_classifier: Union[bool, Classifier] = False,
        scale_adversarial_loss: Union[float, Literal["auto"]] = "auto",
        **loss_kwargs,
    ):
        super().__init__(
            module=module,
            optimizer=optimizer,
            optimizer_creator=optimizer_creator,
            lr=lr,
            weight_decay=weight_decay,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_threshold=lr_threshold,
            lr_scheduler_metric=lr_scheduler_metric,
            lr_min=lr_min,
            **loss_kwargs,
        )
        if adversarial_classifier is True:
            self.n_output_classifier = self.module.n_batch
            self.adversarial_classifier = Classifier(
                n_input=self.module.n_latent,
                n_hidden=32,
                n_labels=self.n_output_classifier,
                n_layers=2,
                logits=True,
            )
        else:
            self.adversarial_classifier = adversarial_classifier
        self.scale_adversarial_loss = scale_adversarial_loss

    def loss_adversarial_classifier(self, z, batch_index, predict_true_class=True):
        """Loss for adversarial classifier."""
        n_classes = self.n_output_classifier
        cls_logits = torch.nn.LogSoftmax(dim=1)(self.adversarial_classifier(z))

        if predict_true_class:
            cls_target = one_hot(batch_index, n_classes)
        else:
            one_hot_batch = one_hot(batch_index, n_classes)
            cls_target = torch.zeros_like(one_hot_batch)
            # place zeroes where true label is
            cls_target.masked_scatter_(
                ~one_hot_batch.bool(), torch.ones_like(one_hot_batch) / (n_classes - 1)
            )

        l_soft = cls_logits * cls_target
        loss = -l_soft.sum(dim=1).mean()

        return loss

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        """Training step for adversarial training."""
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        batch_tensor = batch[REGISTRY_KEYS.BATCH_KEY]
        if optimizer_idx == 0:
            inference_outputs, _, scvi_loss = self.forward(
                batch, loss_kwargs=self.loss_kwargs
            )
            loss = scvi_loss.loss
            # fool classifier if doing adversarial training
            if kappa > 0 and self.adversarial_classifier is not False:
                z = inference_outputs["z"]
                fool_loss = self.loss_adversarial_classifier(z, batch_tensor, False)
                loss += fool_loss * kappa

            self.log("train_loss", loss, on_epoch=True)
            self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")
            return loss

        # train adversarial classifier
        # this condition will not be met unless self.adversarial_classifier is not False
        if optimizer_idx == 1:
            inference_inputs = self.module._get_inference_input(batch)
            outputs = self.module.inference(**inference_inputs)
            z = outputs["z"]
            loss = self.loss_adversarial_classifier(z.detach(), batch_tensor, True)
            loss *= kappa

            return loss

    def configure_optimizers(self):
        """Configure optimizers for adversarial training."""
        params1 = filter(lambda p: p.requires_grad, self.module.parameters())
        optimizer1 = self.get_optimizer_creator()(params1)
        config1 = {"optimizer": optimizer1}
        if self.reduce_lr_on_plateau:
            scheduler1 = ReduceLROnPlateau(
                optimizer1,
                patience=self.lr_patience,
                factor=self.lr_factor,
                threshold=self.lr_threshold,
                min_lr=self.lr_min,
                threshold_mode="abs",
                verbose=True,
            )
            config1.update(
                {
                    "lr_scheduler": scheduler1,
                    "monitor": self.lr_scheduler_metric,
                },
            )

        if self.adversarial_classifier is not False:
            params2 = filter(
                lambda p: p.requires_grad, self.adversarial_classifier.parameters()
            )
            optimizer2 = torch.optim.Adam(
                params2, lr=1e-3, eps=0.01, weight_decay=self.weight_decay
            )
            config2 = {"optimizer": optimizer2}

            # bug in pytorch lightning requires this way to return
            opts = [config1.pop("optimizer"), config2["optimizer"]]
            if "lr_scheduler" in config1:
                config1["scheduler"] = config1.pop("lr_scheduler")
                scheds = [config1]
                return opts, scheds
            else:
                return opts

        return config1


class SemiSupervisedTrainingPlan(TrainingPlan):
    """
    Lightning module task for SemiSupervised Training.

    Parameters
    ----------
    module
        A module instance from class ``BaseModuleClass``.
    classification_ratio
        Weight of the classification_loss in loss function
    lr
        Learning rate used for optimization :class:`~torch.optim.Adam`.
    weight_decay
        Weight decay used in :class:`~torch.optim.Adam`.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
        Only activated when `n_epochs_kl_warmup` is set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from 0 to 1.
        Overrides `n_steps_kl_warmup` when both are not `None`.
    reduce_lr_on_plateau
        Whether to monitor validation loss and reduce learning rate when validation set
        `lr_scheduler_metric` plateaus.
    lr_factor
        Factor to reduce learning rate.
    lr_patience
        Number of epochs with no improvement after which learning rate will be reduced.
    lr_threshold
        Threshold for measuring the new optimum.
    lr_scheduler_metric
        Which metric to track for learning rate reduction.
    **loss_kwargs
        Keyword args to pass to the loss method of the `module`.
        `kl_weight` should not be passed here and is handled automatically.
    """

    def __init__(
        self,
        module: BaseModuleClass,
        *,
        classification_ratio: int = 50,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        **loss_kwargs,
    ):
        super().__init__(
            module=module,
            lr=lr,
            weight_decay=weight_decay,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_threshold=lr_threshold,
            lr_scheduler_metric=lr_scheduler_metric,
            **loss_kwargs,
        )
        self.loss_kwargs.update({"classification_ratio": classification_ratio})

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        """Training step for semi-supervised training."""
        # Potentially dangerous if batch is from a single dataloader with two keys
        if len(batch) == 2:
            full_dataset = batch[0]
            labelled_dataset = batch[1]
        else:
            full_dataset = batch
            labelled_dataset = None

        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        input_kwargs = dict(
            feed_labels=False,
            labelled_tensors=labelled_dataset,
        )
        input_kwargs.update(self.loss_kwargs)
        _, _, loss_output = self.forward(full_dataset, loss_kwargs=input_kwargs)
        loss = loss_output.loss
        self.log(
            "train_loss",
            loss,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.compute_and_log_metrics(loss_output, self.train_metrics, "train")
        return loss

    def validation_step(self, batch, batch_idx, optimizer_idx=0):
        """Validation step for semi-supervised training."""
        # Potentially dangerous if batch is from a single dataloader with two keys
        if len(batch) == 2:
            full_dataset = batch[0]
            labelled_dataset = batch[1]
        else:
            full_dataset = batch
            labelled_dataset = None

        input_kwargs = dict(
            feed_labels=False,
            labelled_tensors=labelled_dataset,
        )
        input_kwargs.update(self.loss_kwargs)
        _, _, loss_output = self.forward(full_dataset, loss_kwargs=input_kwargs)
        loss = loss_output.loss
        self.log(
            "validation_loss",
            loss,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.compute_and_log_metrics(loss_output, self.val_metrics, "validation")


class LowLevelPyroTrainingPlan(TunableMixin, pl.LightningModule):
    """
    Lightning module task to train Pyro scvi-tools modules.

    Parameters
    ----------
    pyro_module
        An instance of :class:`~scvi.module.base.PyroBaseModuleClass`. This object
        should have callable `model` and `guide` attributes or methods.
    loss_fn
        A Pyro loss. Should be a subclass of :class:`~pyro.infer.ELBO`.
        If `None`, defaults to :class:`~pyro.infer.Trace_ELBO`.
    optim
        A Pytorch optimizer class, e.g., :class:`~torch.optim.Adam`. If `None`,
        defaults to :class:`torch.optim.Adam`.
    optim_kwargs
        Keyword arguments for optimiser. If `None`, defaults to `dict(lr=1e-3)`.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
        Only activated when `n_epochs_kl_warmup` is set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from 0 to 1.
        Overrides `n_steps_kl_warmup` when both are not `None`.
    scale_elbo
        Scale ELBO using :class:`~pyro.poutine.scale`. Potentially useful for avoiding
        numerical inaccuracy when working with very large ELBO.
    """

    def __init__(
        self,
        pyro_module: PyroBaseModuleClass,
        loss_fn: Optional[pyro.infer.ELBO] = None,
        optim: Optional[torch.optim.Adam] = None,
        optim_kwargs: Optional[dict] = None,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        scale_elbo: float = 1.0,
    ):
        super().__init__()
        self.module = pyro_module
        self._n_obs_training = None

        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()
        if "lr" not in optim_kwargs.keys():
            optim_kwargs.update({"lr": 1e-3})
        self.optim_kwargs = optim_kwargs

        self.loss_fn = pyro.infer.Trace_ELBO() if loss_fn is None else loss_fn
        self.optim = torch.optim.Adam if optim is None else optim
        self.n_steps_kl_warmup = n_steps_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.use_kl_weight = False
        if isinstance(self.module.model, PyroModule):
            self.use_kl_weight = (
                "kl_weight" in signature(self.module.model.forward).parameters
            )
        elif callable(self.module.model):
            self.use_kl_weight = "kl_weight" in signature(self.module.model).parameters
        self.scale_elbo = scale_elbo
        self.scale_fn = (
            lambda obj: pyro.poutine.scale(obj, self.scale_elbo)
            if self.scale_elbo != 1
            else obj
        )
        self.differentiable_loss_fn = self.loss_fn.differentiable_loss

    def training_step(self, batch, batch_idx):
        """Training step for Pyro training."""
        args, kwargs = self.module._get_fn_args_from_batch(batch)
        # Set KL weight if necessary.
        # Note: if applied, ELBO loss in progress bar is the effective KL annealed loss, not the true ELBO.
        if self.use_kl_weight:
            kwargs.update({"kl_weight": self.kl_weight})
        # pytorch lightning requires a Tensor object for loss
        loss = self.differentiable_loss_fn(
            self.scale_fn(self.module.model),
            self.scale_fn(self.module.guide),
            *args,
            **kwargs,
        )
        return {"loss": loss}

    def training_epoch_end(self, outputs):
        """Training epoch end for Pyro training."""
        elbo = 0
        n = 0
        for out in outputs:
            elbo += out["loss"]
            n += 1
        elbo /= n
        self.log("elbo_train", elbo, prog_bar=True)

    def configure_optimizers(self):
        """Configure optimizers for the model."""
        return self.optim(self.module.parameters(), **self.optim_kwargs)

    def forward(self, *args, **kwargs):
        """Passthrough to the model's forward method."""
        return self.module(*args, **kwargs)

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        return _compute_kl_weight(
            self.current_epoch,
            self.global_step,
            self.n_epochs_kl_warmup,
            self.n_steps_kl_warmup,
            min_kl_weight=1e-3,
        )

    @property
    def n_obs_training(self):
        """
        Number of training examples.

        If not `None`, updates the `n_obs` attr
        of the Pyro module's `model` and `guide`, if they exist.
        """
        return self._n_obs_training

    @n_obs_training.setter
    def n_obs_training(self, n_obs: int):
        # important for scaling log prob in Pyro plates
        if n_obs is not None:
            if hasattr(self.module.model, "n_obs"):
                self.module.model.n_obs = n_obs
            if hasattr(self.module.guide, "n_obs"):
                self.module.guide.n_obs = n_obs

        self._n_obs_training = n_obs


class PyroTrainingPlan(LowLevelPyroTrainingPlan):
    """
    Lightning module task to train Pyro scvi-tools modules.

    Parameters
    ----------
    pyro_module
        An instance of :class:`~scvi.module.base.PyroBaseModuleClass`. This object
        should have callable `model` and `guide` attributes or methods.
    loss_fn
        A Pyro loss. Should be a subclass of :class:`~pyro.infer.ELBO`.
        If `None`, defaults to :class:`~pyro.infer.Trace_ELBO`.
    optim
        A Pyro optimizer instance, e.g., :class:`~pyro.optim.Adam`. If `None`,
        defaults to :class:`pyro.optim.Adam` optimizer with a learning rate of `1e-3`.
    optim_kwargs
        Keyword arguments for **default** optimiser :class:`pyro.optim.Adam`.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
        Only activated when `n_epochs_kl_warmup` is set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from 0 to 1.
        Overrides `n_steps_kl_warmup` when both are not `None`.
    scale_elbo
        Scale ELBO using :class:`~pyro.poutine.scale`. Potentially useful for avoiding
        numerical inaccuracy when working with very large ELBO.
    """

    def __init__(
        self,
        pyro_module: PyroBaseModuleClass,
        loss_fn: Optional[pyro.infer.ELBO] = None,
        optim: Optional[pyro.optim.PyroOptim] = None,
        optim_kwargs: Optional[dict] = None,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        scale_elbo: float = 1.0,
    ):
        super().__init__(
            pyro_module=pyro_module,
            loss_fn=loss_fn,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            n_steps_kl_warmup=n_steps_kl_warmup,
            scale_elbo=scale_elbo,
        )
        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()
        if "lr" not in optim_kwargs.keys():
            optim_kwargs.update({"lr": 1e-3})
        self.optim = (
            pyro.optim.Adam(optim_args=optim_kwargs) if optim is None else optim
        )
        # We let SVI take care of all optimization
        self.automatic_optimization = False

        self.svi = pyro.infer.SVI(
            model=self.scale_fn(self.module.model),
            guide=self.scale_fn(self.module.guide),
            optim=self.optim,
            loss=self.loss_fn,
        )
        # See configure_optimizers for what this does
        self._dummy_param = torch.nn.Parameter(torch.Tensor([0.0]))

    def training_step(self, batch, batch_idx):
        """Training step for Pyro training."""
        args, kwargs = self.module._get_fn_args_from_batch(batch)
        # Set KL weight if necessary.
        # Note: if applied, ELBO loss in progress bar is the effective KL annealed loss, not the true ELBO.
        if self.use_kl_weight:
            kwargs.update({"kl_weight": self.kl_weight})
        # pytorch lightning requires a Tensor object for loss
        loss = torch.Tensor([self.svi.step(*args, **kwargs)])

        _opt = self.optimizers()
        _opt.step()

        return {"loss": loss}

    def configure_optimizers(self):
        """
        Shim optimizer for PyTorch Lightning.

        PyTorch Lightning wants to take steps on an optimizer
        returned by this function in order to increment the global
        step count. See PyTorch Lighinting optimizer manual loop.

        Here we provide a shim optimizer that we can take steps on
        at minimal computational cost in order to keep Lightning happy :).
        """
        return torch.optim.Adam([self._dummy_param])

    def optimizer_step(self, *args, **kwargs):  # noqa: D102
        pass

    def backward(self, *args, **kwargs):  # noqa: D102
        pass


class ClassifierTrainingPlan(TunableMixin, pl.LightningModule):
    """
    Lightning module task to train a simple MLP classifier.

    Parameters
    ----------
    classifier
        A model instance from :class:`~scvi.module.Classifier`.
    lr
        Learning rate used for optimization.
    weight_decay
        Weight decay used in optimization.
    eps
        eps used for optimization.
    optimizer
        One of "Adam" (:class:`~torch.optim.Adam`), "AdamW" (:class:`~torch.optim.AdamW`).
    data_key
        Key for classifier input in tensor dict minibatch
    labels_key
        Key for classifier label in tensor dict minibatch
    loss
        PyTorch loss to use
    """

    def __init__(
        self,
        classifier: BaseModuleClass,
        *,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        eps: float = 0.01,
        optimizer: Literal["Adam", "AdamW"] = "Adam",
        data_key: str = REGISTRY_KEYS.X_KEY,
        labels_key: str = REGISTRY_KEYS.LABELS_KEY,
        loss: Callable = torch.nn.CrossEntropyLoss,
    ):
        super().__init__()
        self.module = classifier
        self.lr = lr
        self.weight_decay = weight_decay
        self.eps = eps
        self.optimizer_name = optimizer
        self.data_key = data_key
        self.labels_key = labels_key
        self.loss_fn = loss()

        if self.module.logits is False and loss == torch.nn.CrossEntropyLoss:
            raise UserWarning(
                "classifier should return logits when using CrossEntropyLoss."
            )

    def forward(self, *args, **kwargs):
        """Passthrough to the module's forward function."""
        return self.module(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        """Training step for classifier training."""
        soft_prediction = self.forward(batch[self.data_key])
        loss = self.loss_fn(soft_prediction, batch[self.labels_key].view(-1).long())
        self.log("train_loss", loss, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        """Validation step for classifier training."""
        soft_prediction = self.forward(batch[self.data_key])
        loss = self.loss_fn(soft_prediction, batch[self.labels_key].view(-1).long())
        self.log("validation_loss", loss)

        return loss

    def configure_optimizers(self):
        """Configure optimizers for classifier training."""
        params = filter(lambda p: p.requires_grad, self.module.parameters())
        if self.optimizer_name == "Adam":
            optim_cls = torch.optim.Adam
        elif self.optimizer_name == "AdamW":
            optim_cls = torch.optim.AdamW
        else:
            raise ValueError("Optimizer not understood.")
        optimizer = optim_cls(
            params, lr=self.lr, eps=self.eps, weight_decay=self.weight_decay
        )

        return optimizer


class JaxTrainingPlan(TrainingPlan):
    """
    Lightning module task to train Pyro scvi-tools modules.

    Parameters
    ----------
    module
        An instance of :class:`~scvi.module.base.JaxModuleWraper`.
    optimizer
        One of "Adam", "AdamW", or "Custom", which requires a custom
        optimizer creator callable to be passed via `optimizer_creator`.
    optimizer_creator
        A callable returning a :class:`~optax.GradientTransformation`.
        This allows using any optax optimizer with custom hyperparameters.
    lr
        Learning rate used for optimization, when `optimizer_creator` is None.
    weight_decay
        Weight decay used in optimization, when `optimizer_creator` is None.
    eps
        eps used for optimization, when `optimizer_creator` is None.
    max_norm
        Max global norm of gradients for gradient clipping.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`. Only activated when `n_epochs_kl_warmup` is
        set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from `min_kl_weight` to
        `max_kl_weight`. Overrides `n_steps_kl_warmup` when both are not `None`.
    """

    def __init__(
        self,
        module: JaxBaseModuleClass,
        *,
        optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam",
        optimizer_creator: Optional[JaxOptimizerCreator] = None,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        eps: float = 0.01,
        max_norm: Optional[float] = None,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        **loss_kwargs,
    ):
        super().__init__(
            module=module,
            lr=lr,
            weight_decay=weight_decay,
            eps=eps,
            optimizer=optimizer,
            optimizer_creator=optimizer_creator,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            **loss_kwargs,
        )
        self.max_norm = max_norm
        self.automatic_optimization = False
        self._dummy_param = torch.nn.Parameter(torch.Tensor([0.0]))

    def get_optimizer_creator(self) -> JaxOptimizerCreator:
        """Get optimizer creator for the model."""
        clip_by = (
            optax.clip_by_global_norm(self.max_norm)
            if self.max_norm
            else optax.identity()
        )
        if self.optimizer_name == "Adam":
            # Replicates PyTorch Adam defaults
            optim = optax.chain(
                clip_by,
                optax.additive_weight_decay(weight_decay=self.weight_decay),
                optax.adam(self.lr, eps=self.eps),
            )
        elif self.optimizer_name == "AdamW":
            optim = optax.chain(
                clip_by,
                optax.clip_by_global_norm(self.max_norm),
                optax.adamw(self.lr, eps=self.eps, weight_decay=self.weight_decay),
            )
        elif self.optimizer_name == "Custom":
            optim = self._optimizer_creator
        else:
            raise ValueError("Optimizer not understood.")

        return lambda: optim

    def set_train_state(self, params, state=None):
        """Set the state of the module."""
        if self.module.train_state is not None:
            return
        optimizer = self.get_optimizer_creator()()
        train_state = TrainStateWithState.create(
            apply_fn=self.module.apply,
            params=params,
            tx=optimizer,
            state=state,
        )
        self.module.train_state = train_state

    @staticmethod
    @jax.jit
    def jit_training_step(
        state: TrainStateWithState,
        batch: Dict[str, np.ndarray],
        rngs: Dict[str, jnp.ndarray],
        **kwargs,
    ):
        """Jit training step."""
        # state can't be passed here
        def loss_fn(params):
            vars_in = {"params": params, **state.state}
            outputs, new_model_state = state.apply_fn(
                vars_in, batch, rngs=rngs, mutable=list(state.state.keys()), **kwargs
            )
            loss_output = outputs[2]
            loss = loss_output.loss
            return loss, (loss_output, new_model_state)

        (loss, (loss_output, new_model_state)), grads = jax.value_and_grad(
            loss_fn, has_aux=True
        )(state.params)
        new_state = state.apply_gradients(grads=grads, state=new_model_state)
        return new_state, loss, loss_output

    def training_step(self, batch, batch_idx):
        """Training step for Jax."""
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        self.module.train()
        self.module.train_state, _, loss_output = self.jit_training_step(
            self.module.train_state,
            batch,
            self.module.rngs,
            loss_kwargs=self.loss_kwargs,
        )
        loss_output = jax.tree_util.tree_map(
            lambda x: torch.tensor(jax.device_get(x)),
            loss_output,
        )
        # TODO: Better way to get batch size
        self.log(
            "train_loss",
            loss_output.loss,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
            prog_bar=True,
        )
        self.compute_and_log_metrics(loss_output, self.train_metrics, "train")
        # Update the dummy optimizer to update the global step
        _opt = self.optimizers()
        _opt.step()

    @partial(jax.jit, static_argnums=(0,))
    def jit_validation_step(
        self,
        state: TrainStateWithState,
        batch: Dict[str, np.ndarray],
        rngs: Dict[str, jnp.ndarray],
        **kwargs,
    ):
        """Jit validation step."""
        vars_in = {"params": state.params, **state.state}
        outputs = self.module.apply(vars_in, batch, rngs=rngs, **kwargs)
        loss_output = outputs[2]

        return loss_output

    def validation_step(self, batch, batch_idx):
        """Validation step for Jax."""
        self.module.eval()
        loss_output = self.jit_validation_step(
            self.module.train_state,
            batch,
            self.module.rngs,
            loss_kwargs=self.loss_kwargs,
        )
        loss_output = jax.tree_util.tree_map(
            lambda x: torch.tensor(jax.device_get(x)),
            loss_output,
        )
        self.log(
            "validation_loss",
            loss_output.loss,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.compute_and_log_metrics(loss_output, self.val_metrics, "validation")

    @staticmethod
    def transfer_batch_to_device(batch, device, dataloader_idx):
        """Bypass Pytorch Lightning device management."""
        return batch

    def configure_optimizers(self):
        """
        Shim optimizer for PyTorch Lightning.

        PyTorch Lightning wants to take steps on an optimizer
        returned by this function in order to increment the global
        step count. See PyTorch Lighinting optimizer manual loop.

        Here we provide a shim optimizer that we can take steps on
        at minimal computational cost in order to keep Lightning happy :).
        """
        return torch.optim.Adam([self._dummy_param])

    def optimizer_step(self, *args, **kwargs):  # noqa: D102
        pass

    def backward(self, *args, **kwargs):  # noqa: D102
        pass

    def forward(self, *args, **kwargs):  # noqa: D102
        pass
