import torch
from pytorch_lightning.loggers import LightningLoggerBase
from pytorch_lightning.utilities import rank_zero_only


class SimpleLogger(LightningLoggerBase):
    def __init__(self):
        super().__init__()
        self._data = {}

    def experiment(self):
        # Return the experiment object associated with this logger.
        pass

    @rank_zero_only
    def log_hyperparams(self, params):
        # params is an argparse.Namespace
        # your code to record hyperparameters goes here
        pass

    @rank_zero_only
    def log_metrics(self, metrics, step):
        """Record metrics."""

        def _handle_value(value):
            if isinstance(value, torch.Tensor):
                return value.item()
            return value

        for metric, value in metrics.items():
            if metric not in self._data:
                self._data[metric] = []
            self._data[metric].append(_handle_value(value))

    @property
    def history(self):
        return self._data

    @property
    def version(self):
        return "1"

    @property
    def name(self):
        return "SimpleLogger"
