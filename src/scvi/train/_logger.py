import os
import pickle
from typing import Any

import pandas as pd
import torch
from lightning.pytorch.loggers.logger import Logger, rank_zero_experiment
from lightning.pytorch.utilities import rank_zero_only


class SimpleExperiment:
    """Simple experiment class."""

    def __init__(self):
        self.data = {}

    def log_hparams(self, params: dict[str, Any]) -> None:
        """Record hparams."""

    def log_metrics(self, metrics: dict[str, float], step: int | None = None) -> None:
        """Record metrics."""

        def _handle_value(value):
            if isinstance(value, torch.Tensor):
                return value.item()
            return value

        if "step" in metrics.keys():
            time_point = metrics.pop("step")
            time_point_name = "step"
        elif "epoch" in metrics.keys():
            time_point = metrics.pop("epoch")
            time_point_name = "epoch"
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
        self,
        name: str = "lightning_logs",
        version: int | str | None = None,
        save_dir: str | None = None,
    ):
        super().__init__()
        self._name = name
        self._experiment = None
        self._version = version
        self._save_dir = save_dir or os.getcwd()
        # run directory like: <save_dir>/<name>/version_<N>
        self._run_dir = os.path.join(self._save_dir, self._name, f"version_{self.version}")
        os.makedirs(self._run_dir, exist_ok=True)
        self.history_path = os.path.join(self._run_dir, "history.pkl")

    @property
    @rank_zero_experiment
    def experiment(self):
        """Return the experiment object associated with this logger."""
        if self._experiment is None:
            self._experiment = SimpleExperiment()
        return self._experiment

    @rank_zero_only
    def log_hyperparams(self, params):
        # params is an argparse.Namespace
        # your code to record hyperparameters goes here
        pass

    @rank_zero_only
    def log_metrics(self, metrics, step):
        """Record metrics."""
        self.experiment.log_metrics(metrics, step)

    @property
    def history(self) -> dict[str, pd.DataFrame]:
        return getattr(self.experiment, "data", {})  # {} instead of AttributeError on non-rank0

    @rank_zero_only
    def finalize(self, status: str) -> None:
        # Persist history from rank-0 AFTER training ends
        try:
            with open(self.history_path, "wb") as f:
                pickle.dump(self.history, f)
        except (OSError, pickle.PickleError) as e:
            print(f"[SimpleLogger] Failed to save history: {e}")

    @property
    def version(self) -> int:
        """Gets the version of the experiment.

        Returns
        -------
        The version of the experiment if it is specified, else the next version.
        """
        if self._version is None:
            self._version = 1
        return self._version

    @property
    def name(self):
        return self._name

    @property
    def save_dir(self):
        return self._save_dir
