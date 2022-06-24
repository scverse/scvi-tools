from functools import partial
from inspect import getfullargspec, signature
from typing import Callable, Dict, Optional, Union

import jax
import jax.numpy as jnp
import numpy as np
import optax
import pyro
import pytorch_lightning as pl
import torch
from pyro.nn import PyroModule
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torchmetrics import MetricCollection

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.module import Classifier
from scvi.module.base import (
    BaseModuleClass,
    JaxModuleWrapper,
    LossRecorder,
    PyroBaseModuleClass,
    TrainStateWithBatchNorm,
)
from scvi.nn import one_hot

from ._metrics import ElboMetric


def _compute_kl_weight(
    epoch: int,
    step: int,
    n_epochs_kl_warmup: Optional[int],
    n_steps_kl_warmup: Optional[int],
    min_weight: Optional[float] = None,
) -> float:
    epoch_criterion = n_epochs_kl_warmup is not None
    step_criterion = n_steps_kl_warmup is not None
    if epoch_criterion:
        kl_weight = min(1.0, epoch / n_epochs_kl_warmup)
    elif step_criterion:
        kl_weight = min(1.0, step / n_steps_kl_warmup)
    else:
        kl_weight = 1.0
    if min_weight is not None:
        kl_weight = max(kl_weight, min_weight)
    return kl_weight


class TrainingPlan(pl.LightningModule):
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
    lr

        Learning rate used for optimization.
    weight_decay
        Weight decay used in optimizatoin.
    eps
        eps used for optimization.
    optimizer
        One of "Adam" (:class:`~torch.optim.Adam`), "AdamW" (:class:`~torch.optim.AdamW`).
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
    **loss_kwargs
        Keyword args to pass to the loss method of the `module`.
        `kl_weight` should not be passed here and is handled automatically.
    """

    def __init__(
        self,
        module: BaseModuleClass,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        eps: float = 0.01,
        optimizer: Literal["Adam", "AdamW"] = "Adam",
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        lr_min: float = 0,
        **loss_kwargs,
    ):
        super(TrainingPlan, self).__init__()
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

        self._n_obs_training = None
        self._n_obs_validation = None

        # automatic handling of kl weight
        self._loss_args = getfullargspec(self.module.loss)[0]
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
        collection = MetricCollection(
            {metric.name: metric for metric in [elbo, rec_loss, kl_local, kl_global]}
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
        """Passthrough to `model.forward()`."""
        return self.module(*args, **kwargs)

    @torch.no_grad()
    def compute_and_log_metrics(
        self,
        loss_recorder: LossRecorder,
        metrics: MetricCollection,
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
        rec_loss = loss_recorder.reconstruction_loss
        n_obs_minibatch = rec_loss.shape[0]
        rec_loss = rec_loss.sum()
        kl_local = loss_recorder.kl_local.sum()
        kl_global = loss_recorder.kl_global

        # use the torchmetric object for the ELBO
        metrics.update(
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
        for extra_metric in loss_recorder.extra_metric_attrs:
            met = getattr(loss_recorder, extra_metric)
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{extra_metric}_{mode}",
                met,
                on_step=False,
                on_epoch=True,
                batch_size=n_obs_minibatch,
            )

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("train_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")
        return scvi_loss.loss

    def validation_step(self, batch, batch_idx):
        # loss kwargs here contains `n_obs` equal to n_training_obs
        # so when relevant, the actual loss value is rescaled to number
        # of training examples
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        self.log("validation_loss", scvi_loss.loss, on_epoch=True)
        self.compute_and_log_metrics(scvi_loss, self.val_metrics, "validation")

    def configure_optimizers(self):
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
        )


class AdversarialTrainingPlan(TrainingPlan):
    """
    Train vaes with adversarial loss option to encourage latent space mixing.

    Parameters
    ----------
    module
        A module instance from class ``BaseModuleClass``.
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
        lr=1e-3,
        weight_decay=1e-6,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
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
        params1 = filter(lambda p: p.requires_grad, self.module.parameters())
        optimizer1 = torch.optim.Adam(
            params1, lr=self.lr, eps=0.01, weight_decay=self.weight_decay
        )
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
        classification_ratio: int = 50,
        lr=1e-3,
        weight_decay=1e-6,
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
        super(SemiSupervisedTrainingPlan, self).__init__(
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
        _, _, scvi_losses = self.forward(full_dataset, loss_kwargs=input_kwargs)
        loss = scvi_losses.loss
        self.log(
            "train_loss",
            loss,
            on_epoch=True,
            batch_size=len(scvi_losses.reconstruction_loss),
        )
        self.compute_and_log_metrics(scvi_losses, self.train_metrics, "train")
        return loss

    def validation_step(self, batch, batch_idx, optimizer_idx=0):
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
        _, _, scvi_losses = self.forward(full_dataset, loss_kwargs=input_kwargs)
        loss = scvi_losses.loss
        self.log(
            "validation_loss",
            loss,
            on_epoch=True,
            batch_size=len(scvi_losses.reconstruction_loss),
        )
        self.compute_and_log_metrics(scvi_losses, self.val_metrics, "validation")


class PyroTrainingPlan(pl.LightningModule):
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
        super().__init__()
        self.module = pyro_module
        self._n_obs_training = None

        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()
        if "lr" not in optim_kwargs.keys():
            optim_kwargs.update({"lr": 1e-3})

        self.loss_fn = pyro.infer.Trace_ELBO() if loss_fn is None else loss_fn
        self.optim = (
            pyro.optim.Adam(optim_args=optim_kwargs) if optim is None else optim
        )

        self.n_steps_kl_warmup = n_steps_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup

        self.automatic_optimization = False

        self.use_kl_weight = False
        if isinstance(self.module.model, PyroModule):
            self.use_kl_weight = (
                "kl_weight" in signature(self.module.model.forward).parameters
            )

        elif callable(self.module.model):
            self.use_kl_weight = "kl_weight" in signature(self.module.model).parameters

        def scale(pyro_obj):
            if scale_elbo == 1:
                return pyro_obj
            else:
                return pyro.poutine.scale(pyro_obj, scale_elbo)

        self.svi = pyro.infer.SVI(
            model=scale(self.module.model),
            guide=scale(self.module.guide),
            optim=self.optim,
            loss=self.loss_fn,
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
                setattr(self.module.model, "n_obs", n_obs)
            if hasattr(self.module.guide, "n_obs"):
                setattr(self.module.guide, "n_obs", n_obs)

        self._n_obs_training = n_obs

    def forward(self, *args, **kwargs):
        """Passthrough to `model.forward()`."""
        return self.module(*args, **kwargs)

    def training_step(self, batch, batch_idx):
        args, kwargs = self.module._get_fn_args_from_batch(batch)
        # Set KL weight if necessary.
        # Note: if applied, ELBO loss in progress bar is the effective KL annealed loss, not the true ELBO.
        if self.use_kl_weight:
            kwargs.update({"kl_weight": self.kl_weight})
        # pytorch lightning requires a Tensor object for loss
        loss = torch.Tensor([self.svi.step(*args, **kwargs)])

        return {"loss": loss}

    def training_epoch_end(self, outputs):
        elbo = 0
        n = 0
        for out in outputs:
            elbo += out["loss"]
            n += 1
        elbo /= n
        self.log("elbo_train", elbo, prog_bar=True)

    def configure_optimizers(self):
        return None

    def optimizer_step(self, *args, **kwargs):
        pass

    def backward(self, *args, **kwargs):
        pass

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        return _compute_kl_weight(
            self.current_epoch,
            self.global_step,
            self.n_epochs_kl_warmup,
            self.n_steps_kl_warmup,
            min_weight=1e-3,
        )


class ClassifierTrainingPlan(pl.LightningModule):
    """
    Lightning module task to train a simple MLP classifier.

    Parameters
    ----------
    classifier
        A model instance from :class:`~scvi.module.Classifier`.
    lr
        Learning rate used for optimization.
    weight_decay
        Weight decay used in optimizatoin.
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
        """Passthrough to `model.forward()`."""
        return self.module(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):

        soft_prediction = self.forward(batch[self.data_key])
        loss = self.loss_fn(soft_prediction, batch[self.labels_key].view(-1).long())
        self.log("train_loss", loss, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        soft_prediction = self.forward(batch[self.data_key])
        loss = self.loss_fn(soft_prediction, batch[self.labels_key].view(-1).long())
        self.log("validation_loss", loss)

        return loss

    def configure_optimizers(self):
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


class JaxTrainingPlan(pl.LightningModule):
    """
    Lightning module task to train Pyro scvi-tools modules.

    Parameters
    ----------
    module
        An instance of :class:`~scvi.module.base.JaxModuleWraper`.
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
        Only activated when `n_epochs_kl_warmup` is set to None.
    n_epochs_kl_warmup
        Number of epochs to scale weight on KL divergences from 0 to 1.
        Overrides `n_steps_kl_warmup` when both are not `None`.
    """

    def __init__(
        self,
        module: JaxModuleWrapper,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        optim_kwargs: Optional[dict] = None,
        **loss_kwargs,
    ):
        super().__init__()
        self.module = module
        self._n_obs_training = None
        self.loss_kwargs = loss_kwargs
        self.n_steps_kl_warmup = n_steps_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup

        self.automatic_optimization = False

        # automatic handling of kl weight
        self._loss_args = signature(self.module.loss).parameters
        if "kl_weight" in self._loss_args:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})

        # set optim kwargs
        self.optim_kwargs = dict(learning_rate=1e-3, eps=0.01, weight_decay=1e-6)
        if optim_kwargs is not None:
            self.optim_kwargs.update(optim_kwargs)

    def set_train_state(self, params, batch_stats=None):
        self.module.train()

        if self.module.train_state is not None:
            return

        weight_decay = self.optim_kwargs.pop("weight_decay")
        # replicates PyTorch Adam
        optimizer = optax.chain(
            optax.additive_weight_decay(weight_decay=weight_decay),
            optax.adam(**self.optim_kwargs),
        )
        train_state = TrainStateWithBatchNorm.create(
            apply_fn=self.module.apply,
            params=params,
            tx=optimizer,
            batch_stats=batch_stats,
        )
        self.module.train_state = train_state

    @staticmethod
    @jax.jit
    def jit_training_step(
        state: TrainStateWithBatchNorm,
        batch: Dict[str, np.ndarray],
        rngs: Dict[str, jnp.ndarray],
        **kwargs,
    ):
        # batch stats can't be passed here
        def loss_fn(params):
            vars_in = {"params": params, "batch_stats": state.batch_stats}
            outputs, new_model_state = state.apply_fn(
                vars_in, batch, rngs=rngs, mutable=["batch_stats"], **kwargs
            )
            loss_recorder = outputs[2]
            loss = loss_recorder.loss
            elbo = jnp.mean(loss_recorder.reconstruction_loss + loss_recorder.kl_local)
            return loss, (elbo, new_model_state)

        (loss, (elbo, new_model_state)), grads = jax.value_and_grad(
            loss_fn, has_aux=True
        )(state.params)
        new_state = state.apply_gradients(
            grads=grads, batch_stats=new_model_state["batch_stats"]
        )
        return new_state, loss, elbo

    def training_step(self, batch, batch_idx):
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        self.module.train()
        self.module.train_state, loss, elbo = self.jit_training_step(
            self.module.train_state,
            batch,
            self.module.rngs,
            loss_kwargs=self.loss_kwargs,
        )
        loss = torch.tensor(jax.device_get(loss))
        elbo = torch.tensor(jax.device_get(elbo))
        # TODO: Better way to get batch size
        self.log(
            "train_loss",
            loss,
            on_epoch=True,
            batch_size=batch[REGISTRY_KEYS.X_KEY].shape[0],
        )
        self.log(
            "elbo_train",
            elbo,
            on_epoch=True,
            batch_size=batch[REGISTRY_KEYS.X_KEY].shape[0],
        )

    @partial(jax.jit, static_argnums=(0,))
    def jit_validation_step(
        self,
        state: TrainStateWithBatchNorm,
        batch: Dict[str, np.ndarray],
        rngs: Dict[str, jnp.ndarray],
        **kwargs,
    ):
        vars_in = {"params": state.params, "batch_stats": state.batch_stats}
        outputs = self.module.apply(vars_in, batch, rngs=rngs, **kwargs)
        loss_recorder = outputs[2]
        loss = loss_recorder.loss
        elbo = jnp.mean(loss_recorder.reconstruction_loss + loss_recorder.kl_local)

        return loss, elbo

    def validation_step(self, batch, batch_idx):
        self.module.eval()
        loss, elbo = self.jit_validation_step(
            self.module.train_state,
            batch,
            self.module.rngs,
            loss_kwargs=self.loss_kwargs,
        )
        loss = torch.tensor(jax.device_get(loss))
        elbo = torch.tensor(jax.device_get(elbo))
        self.log(
            "validation_loss",
            loss,
            on_epoch=True,
            batch_size=batch[REGISTRY_KEYS.X_KEY].shape[0],
        )
        self.log(
            "elbo_validation",
            elbo,
            on_epoch=True,
            batch_size=batch[REGISTRY_KEYS.X_KEY].shape[0],
        )

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        return _compute_kl_weight(
            self.current_epoch,
            self.global_step,
            self.n_epochs_kl_warmup,
            self.n_steps_kl_warmup,
        )

    @staticmethod
    def transfer_batch_to_device(batch, device, dataloader_idx):
        """Bypass Pytorch Lightning device management."""
        return batch

    def configure_optimizers(self):
        return None

    def optimizer_step(self, *args, **kwargs):
        pass

    def backward(self, *args, **kwargs):
        pass

    def forward(self, *args, **kwargs):
        pass
