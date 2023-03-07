from typing import Callable, Dict

import pytorch_lightning as pl

try:
    from ray.tune.integration.pytorch_lightning import TuneCallback
except ImportError:
    pass

from scvi.model.base import BaseModelClass


class TuneMetricCallback(TuneCallback):
    """Callback to report arbitrary metrics to Ray Tune."""

    def __init__(self, model: BaseModelClass, metrics: Dict[str, Callable]):
        super().__init__()
        self._metrics = metrics
        self._model = model

    @property
    def metrics(self):
        """Metrics."""
        return self._metrics

    @property
    def model(self):
        """Model."""
        return self._model

    def _handle(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        # TODO: this is still a work in progress
        if trainer.sanity_checking:
            return

        self.model.is_trained = True

        _metric_vals = {}
        for metric, metric_fn in self.metrics.items():
            _metric_vals[metric] = metric_fn(self.model)

        self.model.is_trained = False
