import torch
from torchmetrics import Metric

from scvi._compat import Literal


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
    dist_sync_on_step
        Synchronize metric state across processes at each ``forward()``
        before returning the value at the step.
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    # Needs to be explicitly set to avoid TorchMetrics UserWarning.
    full_state_update = True
    _N_OBS_MINIBATCH_KEY = "n_obs_minibatch"

    def __init__(
        self,
        name: str,
        mode: Literal["train", "validation"],
        interval: Literal["obs", "batch"],
        dist_sync_on_step: bool = False,
        **kwargs,
    ):
        super().__init__(dist_sync_on_step=dist_sync_on_step, **kwargs)

        self._name = name
        self._mode = mode
        self._interval = interval

        default_val = torch.tensor(0.0)
        self.add_state("elbo_component", default=default_val)
        self.add_state("n_obs", default=default_val)
        self.add_state("n_batches", default=default_val)

    @property
    def mode(self):
        return self._mode

    @property
    def name(self):
        return f"{self._name}_{self.mode}"

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def interval(self):
        return self._interval

    def get_intervals_recorded(self):
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
        if self._N_OBS_MINIBATCH_KEY not in kwargs:
            raise ValueError(
                f"Missing {self._N_OBS_MINIBATCH_KEY} value in metrics update."
            )
        if self._name not in kwargs:
            raise ValueError(f"Missing {self._name} value in metrics update.")

        elbo_component = kwargs[self._name].detach()
        self.elbo_component += elbo_component

        n_obs_minibatch = kwargs[self._N_OBS_MINIBATCH_KEY]
        self.n_obs += n_obs_minibatch
        self.n_batches += 1

    def compute(self):
        return self.elbo_component / self.get_intervals_recorded()
