from pytorch_lightning.callbacks import Callback
from ray import tune
from ray.tune.integration.pytorch_lightning import TuneCallback


class ModelSave(Callback):
    def __init__(self, model):
        super()
        self.model = model

    def on_validation_epoch_end(self, trainer, pl_module, outputs=None):
        if trainer.sanity_checking:
            return
        step = f"epoch={trainer.current_epoch}-step={trainer.global_step}"
        with tune.checkpoint_dir(step=step) as checkpoint_dir:
<<<<<<< HEAD
            self.model.save(checkpoint_dir + "/checkpoint")
=======
            self.model.save(dir_path=checkpoint_dir + "/checkpoint")
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)


class _TuneReportMetricFunctionsCallback(TuneCallback):
    def __init__(
        self,
        metrics=None,
        metric_functions=None,
        on="validation_end",
        model=None,
    ):
        super(_TuneReportMetricFunctionsCallback, self).__init__(on)
        if isinstance(metrics, str):
            metrics = [metrics]
        self._metrics = metrics
        self._metric_functions = metric_functions
        self._model = model

    def _handle(self, trainer, pl_module):
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
