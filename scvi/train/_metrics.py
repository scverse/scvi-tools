import torch
from torchmetrics import Metric

from scvi._compat import Literal

from ._constants import METRIC_KEYS


class ElboMetric(Metric):
    """
    Elbo metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    name
        Name of metric, used as the prefix of the logged name.
    mode
        Train or validation, used as the suffix of the logged name.
    interval
        The interval over which the metric is computed. If "obs", the metric value
        per observation is computed. If "batch", the metric value per batch is computed.
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    # Needs to be explicitly set to avoid TorchMetrics UserWarning.
    full_state_update = False

    def __init__(
        self,
        name: str,
        mode: Literal["train", "validation"],
        interval: Literal["obs", "batch"],
        **kwargs,
    ):
        super().__init__(**kwargs)

        self._name = name
        self._mode = mode
        self._interval = interval

        self.add_state(
            "elbo_component", default=torch.tensor(0.0), dist_reduce_fx="sum"
        )
        self.add_state("n_obs", default=torch.tensor(0.0), dist_reduce_fx="sum")
        self.add_state("n_batches", default=torch.tensor(0.0), dist_reduce_fx="sum")

    @property
    def mode(self):  # noqa: D102
        return self._mode

    @property
    def name(self):  # noqa: D102
        return f"{self._name}_{self.mode}"

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def interval(self):  # noqa: D102
        return self._interval

    def get_intervals_recorded(self):
        """Get intervals recorded."""
        if self.interval == "obs":
            return self.n_obs
        elif self.interval == "batch":
            return self.n_batches
        raise ValueError(f"Unrecognized interval: {self.interval}.")

    def update(
        self,
        **kwargs,
    ):
        """
        Updates this metric for one minibatch.

        Takes kwargs associated with all metrics being updated for a given minibatch.
        Filters for the relevant metric's value and updates this metric.
        """
        if self._name not in kwargs:
            raise ValueError(f"Missing {self._name} value in metrics update.")

        elbo_component = kwargs[self._name]
        self.elbo_component += elbo_component

        n_obs_minibatch = kwargs[METRIC_KEYS.N_OBS_MINIBATCH]
        self.n_obs += n_obs_minibatch
        self.n_batches += 1

    def compute(self):
        """Compute the metric value."""
        return self.elbo_component / self.get_intervals_recorded()
