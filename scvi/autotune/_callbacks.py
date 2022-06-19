import os
from typing import Any, List, Optional, Union

from pytorch_lightning import LightningModule, Trainer
from pytorch_lightning.callbacks import Callback
from ray import tune
from ray.tune.integration.pytorch_lightning import TuneCallback

from scvi.model.base import BaseModelClass


class ModelSave(Callback):
    def __init__(self, model: BaseModelClass):
        super()
        self.model = model

    def on_validation_epoch_end(
        self,
        trainer: Trainer,
        pl_module: Optional[LightningModule],
        outputs: List[Any] = None,
    ):
        if trainer.sanity_checking:
            return
        step = f"epoch={trainer.current_epoch}-step={trainer.global_step}"
        with tune.checkpoint_dir(step=step) as checkpoint_dir:
            self.model.save(dir_path=os.path.join(checkpoint_dir, "checkpoint"))


class _TuneReportMetricFunctionsCallback(TuneCallback):
    def __init__(
        self,
        metrics: str,
        metric_functions: str,
        on: Union[str, List[str]] = "validation_end",
        model: BaseModelClass = None,
    ):
        super(_TuneReportMetricFunctionsCallback, self).__init__(on)
        if isinstance(metrics, str):
            metrics = [metrics]
        self._metrics = metrics
        self._metric_functions = metric_functions
        self._model = model

    def _handle(self, trainer: Trainer, pl_module: Optional[LightningModule]):
        # Don't report if just doing initial validation sanity checks.
        if trainer.sanity_checking:
            return
        if not self._metrics:
            report_dict = {k: v.item() for k, v in trainer.callback_metrics.items()}
        else:
            report_dict = {}
            for key in self._metrics:
                if isinstance(self._metrics, dict):
                    metric = self._metrics[key]
                else:
                    metric = key
                report_dict[key] = trainer.callback_metrics[metric].item()
        if self._metric_functions:
            for key in self._metric_functions:
                report_dict[key] = self._metric_functions[key](self._model)
        tune.report(**report_dict)
