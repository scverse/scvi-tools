from inspect import getfullargspec
from typing import Literal, Union

import torch
from torchmetrics import MetricCollection

from scvi.external.csi.module import LossRecorder
from scvi.module.base import BaseModuleClass
from scvi.train import TrainingPlan

# TODO could make new metric class to not be called elbo metric as used for other metrics as well
from scvi.train._metrics import ElboMetric


class WeightScaling:
    def __init__(
        self,
        weight_start: float,
        weight_end: float,
        point_start: int,
        point_end: int,
        update_on: Literal["epoch", "step"] = "step",
    ):
        """Linearly scale loss weights between start and end weight accordingly to the current training stage

        Parameters
        ----------
        weight_start
            Starting weight value
        weight_end
            End weight vlue
        point_start
            Training point to start scaling - before weight is weight_start
            Since the epochs are counted after they are run,
            the start point must be set to 0 to represent 1st epoch
        point_end
            Training point to end scaling - after weight is weight_end
            Since the epochs are counted after they are run,
            the start point must be set to n-1 to represent the last epoch
        update_on
            Define training progression based on epochs or steps

        """
        self.weight_start = weight_start
        self.weight_end = weight_end
        self.point_start = point_start
        self.point_end = point_end
        if update_on not in ["step", "epoch"]:
            raise ValueError("update_on not recognized")
        self.update_on = update_on

        weight_diff = self.weight_end - self.weight_start
        n_points = self.point_end - self.point_start
        self.slope = weight_diff / n_points

        if (
            self.weight(epoch=self.point_start, step=self.point_start) < 0
            or self.weight(epoch=self.point_end, step=self.point_end) < 0
        ):
            raise ValueError("Specified weight scaling would lead to negative weights")

    def weight(
        self,
        epoch: int,
        step: int,
    ) -> float:
        """
        Computes the weight for the current step/epoch depending on which update type was set in init

        Parameters
        ----------
        epoch
            Current epoch.
        step
            Current step.
        """
        if self.update_on == "epoch":
            point = epoch
        elif self.update_on == "step":
            point = step
        else:
            # This is ensured not to happen by above init check
            raise ValueError("self.update_on not recognised")

        if point < self.point_start:
            return self.weight_start
        elif point > self.point_end:
            return self.weight_end
        else:
            return self.slope * (point - self.point_start) + self.weight_start


class TrainingPlanCustom(TrainingPlan):
    def __init__(
        self,
        module: BaseModuleClass,
        loss_weights: Union[None, dict[str, Union[float, WeightScaling]]] = None,
        log_on_epoch: bool = True,
        log_on_step: bool = False,
        **kwargs,
    ):
        """Extends scvi TrainingPlan for custom support for other losses.

        Parameters
        ----------
        args
            Passed to parent
        log_on_epoch
            See on_epoch of lightning Module log method
        log_on_step
            See on_step of lightning Module log method
        loss_weights
            Specifies how losses should be weighted and how it may change during training
            Dict with keys being loss names and values being loss weights.
            Loss weights can be floats for constant weight or dict of params passed to WeightScaling object
            Note that other loss weight params from the parent class are ignored
            (e.g. n_steps/epochs_kl_warmup and min/max_kl_weight)
        kwargs
            Passed to parent.
            As described in param loss_weights the loss weighting params of parent are ignored
        """
        super().__init__(module, **kwargs)

        self.log_on_epoch = log_on_epoch
        self.log_on_step = log_on_step

        # automatic handling of loss component weights
        if loss_weights is None:
            loss_weights = {}
        # Make weighting object
        for loss, weight in loss_weights.items():
            if isinstance(weight, dict):
                loss_weights[loss] = WeightScaling(**weight)
        self.loss_weights = loss_weights

        # Ensure that all passed loss weight specifications are in available loss params
        # Also update loss kwargs based on specified weights
        self._loss_args = getfullargspec(self.module.loss)[0]
        # Make sure no loss weights are already in loss kwargs (e.g. from parent init)
        for loss in self._loss_args:
            if loss in self.loss_kwargs:
                del self.loss_kwargs[loss]
        for loss, weight in loss_weights.items():
            if loss not in self._loss_args:
                raise ValueError(
                    f"Loss {loss} for which a weight was specified is not in loss parameters"
                )
            # This will also overwrite the kl_weight from parent
            self.loss_kwargs.update({loss: self.compute_loss_weight(weight=weight)})

    def compute_loss_weight(self, weight):
        if isinstance(weight, float):
            return weight
        elif isinstance(weight, int):
            return float(weight)
        elif isinstance(weight, WeightScaling):
            return weight.weight(epoch=self.current_epoch, step=self.global_step)

    @staticmethod
    def _create_elbo_metric_components(
        mode: Literal["train", "validation"], **kwargs
    ) -> (ElboMetric, MetricCollection):
        """
        Initialize the combined loss collection.

        Parameters
        ----------
        mode
            train/validation

        Returns
        -------
        tuple
            Objects for storing the combined loss

        """
        loss = ElboMetric("loss", mode, "obs")
        collection = MetricCollection({metric.name: metric for metric in [loss]})
        return loss, collection

    def initialize_train_metrics(self):
        """Initialize train combined loss.

        TODO could add other losses
        """
        (
            self.loss_train,
            self.train_metrics,
        ) = self._create_elbo_metric_components(
            mode="train", n_total=self.n_obs_training
        )
        self.loss_train.reset()

    def initialize_val_metrics(self):
        """Initialize val combined loss.

        TODO could add other losses
        """
        (
            self.loss_val,
            self.val_metrics,
        ) = self._create_elbo_metric_components(
            mode="validation", n_total=self.n_obs_validation
        )
        self.loss_val.reset()

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
        metrics
            The loss Metric Collection to update
        mode
            Postfix string to add to the metric name of
            extra metrics. If train also logs the loss in progress bar
        """
        n_obs_minibatch = loss_recorder.n_obs
        loss_sum = loss_recorder.loss_sum

        # use the torchmetric object
        metrics.update(
            loss=loss_sum,
            n_obs_minibatch=n_obs_minibatch,
        )

        self.log(
            f"loss_{mode}",
            loss_recorder.loss_sum,
            on_step=self.log_on_step,
            on_epoch=self.log_on_epoch,
            batch_size=n_obs_minibatch,
            prog_bar=True if mode == "train" else False,
            sync_dist=self.use_sync_dist,
        )

        # accumulate extra metrics passed to loss recorder
        for extra_metric in loss_recorder.extra_metric_attrs:
            met = getattr(loss_recorder, extra_metric)
            if isinstance(met, torch.Tensor):
                if met.shape != torch.Size([]):
                    raise ValueError("Extra tracked metrics should be 0-d tensors.")
                met = met.detach()
            self.log(
                f"{extra_metric}_{mode}",
                met,
                on_step=self.log_on_step,
                on_epoch=self.log_on_epoch,
                batch_size=n_obs_minibatch,
                sync_dist=self.use_sync_dist,
            )

    def training_step(self, batch, batch_idx):
        for loss, weight in self.loss_weights.items():
            self.loss_kwargs.update({loss: self.compute_loss_weight(weight=weight)})
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        # combined loss is logged via compute_and_log_metrics
        self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")
        return scvi_loss.loss

    def validation_step(self, batch, batch_idx):
        _, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        # Combined loss is logged via compute_and_log_metrics
        self.compute_and_log_metrics(scvi_loss, self.val_metrics, "validation")

    @property
    def kl_weight(self):
        # Can not raise not implemented error as used in parent init
        pass
