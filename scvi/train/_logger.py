from typing import Any, Dict, Optional, Union

import pandas as pd
import torch
from pytorch_lightning.loggers.logger import Logger, rank_zero_experiment
from pytorch_lightning.utilities import rank_zero_only


class SimpleExperiment:
    """Simple experiment class."""

    def __init__(self):
        self.data = {}

    def log_hparams(self, params: Dict[str, Any]) -> None:
        """Record hparams."""

    def log_metrics(
        self, metrics: Dict[str, float], step: Optional[int] = None
    ) -> None:
        """Record metrics."""

        def _handle_value(value):
            if isinstance(value, torch.Tensor):
                return value.item()
            return value

        if "epoch" in metrics.keys():
            time_point = metrics.pop("epoch")
            time_point_name = "epoch"
        elif "step" in metrics.keys():
            time_point = metrics.pop("step")
            time_point_name = "step"
        else:
            time_point = step
            time_point_name = "step"
        for metric, value in metrics.items():
            if metric not in self.data:
                self.data[metric] = pd.DataFrame(columns=[metric])
                self.data[metric].index.name = time_point_name
            self.data[metric].loc[time_point, metric] = _handle_value(value)

    def save(self) -> None:
        """Save data."""


class SimpleLogger(Logger):
    """Simple logger class."""

    def __init__(
        self, name: str = "lightning_logs", version: Optional[Union[int, str]] = None
    ):
        super().__init__()
        self._name = name
        self._experiment = None
        self._version = version

    @property
    @rank_zero_experiment
    def experiment(self):
        """Return the experiment object associated with this logger."""
        if self._experiment is None:
            self._experiment = SimpleExperiment()
        return self._experiment

    @rank_zero_only
    def log_hyperparams(self, params):  # noqa: D102
        # params is an argparse.Namespace
        # your code to record hyperparameters goes here
        pass

    @rank_zero_only
    def log_metrics(self, metrics, step):
        """Record metrics."""
        self.experiment.log_metrics(metrics, step)

    @property
    def history(self) -> Dict[str, pd.DataFrame]:  # noqa: D102
        return self.experiment.data

    @property
    def version(self) -> int:
        """
        Gets the version of the experiment.

        Returns
        -------
        The version of the experiment if it is specified, else the next version.
        """
        if self._version is None:
            self._version = 1
        return self._version

    @property
    def name(self):  # noqa: D102
        return self._name
