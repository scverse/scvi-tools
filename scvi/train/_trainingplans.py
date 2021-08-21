from inspect import getfullargspec
from typing import Callable, Optional, Union

import pyro
import pytorch_lightning as pl
import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.module import Classifier
from scvi.module.base import BaseModuleClass, PyroBaseModuleClass
from scvi.nn import one_hot


class TrainingPlan(pl.LightningModule):
    """
    Lightning module task to train scvi-tools modules.

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

        # automatic handling of kl weight
        self._loss_args = getfullargspec(self.module.loss)[0]
        if "kl_weight" in self._loss_args:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})

    @property
    def n_obs_training(self):
        """
        Number of observations in the training set.

        This will update the loss kwargs for loss rescaling.
        """
        return self._n_obs_training

    @n_obs_training.setter
    def n_obs_training(self, n_obs: int):
        if "n_obs" in self._loss_args:
            self.loss_kwargs.update({"n_obs": n_obs})
        self._n_obs_training = n_obs

    def forward(self, *args, **kwargs):
        """Passthrough to `model.forward()`."""
        return self.module(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        reconstruction_loss = scvi_loss.reconstruction_loss
        # pytorch lightning automatically backprops on "loss"
        self.log("train_loss", scvi_loss.loss, on_epoch=True)
        return {
            "loss": scvi_loss.loss,
            "reconstruction_loss_sum": reconstruction_loss.sum(),
            "kl_local_sum": scvi_loss.kl_local.sum(),
            "kl_global": scvi_loss.kl_global,
            "n_obs": reconstruction_loss.shape[0],
        }

    def training_epoch_end(self, outputs):
        n_obs, elbo, rec_loss, kl_local = 0, 0, 0, 0
        for tensors in outputs:
            elbo += tensors["reconstruction_loss_sum"] + tensors["kl_local_sum"]
            rec_loss += tensors["reconstruction_loss_sum"]
            kl_local += tensors["kl_local_sum"]
            n_obs += tensors["n_obs"]
        # kl global same for each minibatch
        kl_global = outputs[0]["kl_global"]
        elbo += kl_global
        self.log("elbo_train", elbo / n_obs)
        self.log("reconstruction_loss_train", rec_loss / n_obs)
        self.log("kl_local_train", kl_local / n_obs)
        self.log("kl_global_train", kl_global)

    def validation_step(self, batch, batch_idx):
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        reconstruction_loss = scvi_loss.reconstruction_loss
        self.log("validation_loss", scvi_loss.loss, on_epoch=True)
        return {
            "reconstruction_loss_sum": reconstruction_loss.sum(),
            "kl_local_sum": scvi_loss.kl_local.sum(),
            "kl_global": scvi_loss.kl_global,
            "n_obs": reconstruction_loss.shape[0],
        }

    def validation_epoch_end(self, outputs):
        """Aggregate validation step information."""
        n_obs, elbo, rec_loss, kl_local = 0, 0, 0, 0
        for tensors in outputs:
            elbo += tensors["reconstruction_loss_sum"] + tensors["kl_local_sum"]
            rec_loss += tensors["reconstruction_loss_sum"]
            kl_local += tensors["kl_local_sum"]
            n_obs += tensors["n_obs"]
        # kl global same for each minibatch
        kl_global = outputs[0]["kl_global"]
        elbo += kl_global
        self.log("elbo_validation", elbo / n_obs)
        self.log("reconstruction_loss_validation", rec_loss / n_obs)
        self.log("kl_local_validation", kl_local / n_obs)
        self.log("kl_global_validation", kl_global)

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
        epoch_criterion = self.n_epochs_kl_warmup is not None
        step_criterion = self.n_steps_kl_warmup is not None
        if epoch_criterion:
            kl_weight = min(1.0, self.current_epoch / self.n_epochs_kl_warmup)
        elif step_criterion:
            kl_weight = min(1.0, self.global_step / self.n_steps_kl_warmup)
        else:
            kl_weight = 1.0
        return kl_weight


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
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        batch_tensor = batch[_CONSTANTS.BATCH_KEY]
        if optimizer_idx == 0:
            loss_kwargs = dict(kl_weight=self.kl_weight)
            inference_outputs, _, scvi_loss = self.forward(
                batch, loss_kwargs=loss_kwargs
            )
            loss = scvi_loss.loss
            # fool classifier if doing adversarial training
            if kappa > 0 and self.adversarial_classifier is not False:
                z = inference_outputs["z"]
                fool_loss = self.loss_adversarial_classifier(z, batch_tensor, False)
                loss += fool_loss * kappa

            reconstruction_loss = scvi_loss.reconstruction_loss
            self.log("train_loss", loss, on_epoch=True)
            return {
                "loss": loss,
                "reconstruction_loss_sum": reconstruction_loss.sum(),
                "kl_local_sum": scvi_loss.kl_local.sum(),
                "kl_global": scvi_loss.kl_global,
                "n_obs": reconstruction_loss.shape[0],
            }

        # train adversarial classifier
        # this condition will not be met unless self.adversarial_classifier is not False
        if optimizer_idx == 1:
            inference_inputs = self.module._get_inference_input(batch)
            outputs = self.module.inference(**inference_inputs)
            z = outputs["z"]
            loss = self.loss_adversarial_classifier(z.detach(), batch_tensor, True)
            loss *= kappa

            return loss

    def training_epoch_end(self, outputs):
        # only report from optimizer one loss signature
        if self.adversarial_classifier:
            super().training_epoch_end(outputs[0])
        else:
            super().training_epoch_end(outputs)

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
        reconstruction_loss = scvi_losses.reconstruction_loss
        self.log("train_loss", loss, on_epoch=True)
        loss_dict = {
            "loss": loss,
            "reconstruction_loss_sum": reconstruction_loss.sum(),
            "kl_local_sum": scvi_losses.kl_local.sum(),
            "kl_global": scvi_losses.kl_global,
            "n_obs": reconstruction_loss.shape[0],
        }
        if hasattr(scvi_losses, "classification_loss"):
            loss_dict["classification_loss"] = scvi_losses.classification_loss
            loss_dict["n_labelled_tensors"] = scvi_losses.n_labelled_tensors
        return loss_dict

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
        reconstruction_loss = scvi_losses.reconstruction_loss
        self.log("validation_loss", loss, on_epoch=True)
        loss_dict = {
            "loss": loss,
            "reconstruction_loss_sum": reconstruction_loss.sum(),
            "kl_local_sum": scvi_losses.kl_local.sum(),
            "kl_global": scvi_losses.kl_global,
            "n_obs": reconstruction_loss.shape[0],
        }
        if hasattr(scvi_losses, "classification_loss"):
            loss_dict["classification_loss"] = scvi_losses.classification_loss
            loss_dict["n_labelled_tensors"] = scvi_losses.n_labelled_tensors
        return loss_dict

    def training_epoch_end(self, outputs):
        super().training_epoch_end(outputs)
        classifier_loss, total_labelled_tensors = 0, 0

        for tensors in outputs:
            if "classification_loss" in tensors.keys():
                n_labelled = tensors["n_labelled_tensors"]
                total_labelled_tensors += n_labelled
                classification_loss = tensors["classification_loss"]
                classifier_loss += classification_loss * n_labelled

        if total_labelled_tensors > 0:
            self.log(
                "classification_loss_train", classifier_loss / total_labelled_tensors
            )

    def validation_epoch_end(self, outputs):
        super().validation_epoch_end(outputs)
        classifier_loss, total_labelled_tensors = 0, 0

        for tensors in outputs:
            if "classification_loss" in tensors.keys():
                n_labelled = tensors["n_labelled_tensors"]
                total_labelled_tensors += n_labelled
                classification_loss = tensors["classification_loss"]
                classifier_loss += classification_loss * n_labelled

        if total_labelled_tensors > 0:
            self.log(
                "classification_loss_validation",
                classifier_loss / total_labelled_tensors,
            )


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
    """

    def __init__(
        self,
        pyro_module: PyroBaseModuleClass,
        loss_fn: Optional[pyro.infer.ELBO] = None,
        optim: Optional[pyro.optim.PyroOptim] = None,
        optim_kwargs: Optional[dict] = None,
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

        self.automatic_optimization = False
        self.pyro_guide = self.module.guide
        self.pyro_model = self.module.model

        self.svi = pyro.infer.SVI(
            model=self.pyro_model,
            guide=self.pyro_guide,
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
        loss = self.svi.step(*args, **kwargs)

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
        data_key: str = _CONSTANTS.X_KEY,
        labels_key: str = _CONSTANTS.LABELS_KEY,
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
