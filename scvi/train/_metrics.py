import torch
from torchmetrics import Metric

from scvi._compat import Literal


class BaseElboMetric(Metric):
    """
    Elbo metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    mode
        Train or validation, used for logging names
    dist_sync_on_step
        optional, by default False
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    def __init__(
        self,
        mode: Literal["train", "validation"],
        dist_sync_on_step: bool = False,
        **kwargs
    ):
        super().__init__(dist_sync_on_step=dist_sync_on_step, **kwargs)

        self._mode = mode

        default_val = torch.tensor(0.0)
        self.add_state("elbo_component", default=default_val)
        self.add_state("n_obs", default=default_val)
        self.add_state("n_batches", default=default_val)

    @property
    def mode(self):
        return self._mode

    def update(
        self,
        elbo_component: torch.Tensor,
        n_obs_minibatch: int,
    ):
        """Updates all metrics."""
        elbo_component = elbo_component.detach()
        self.elbo_component += elbo_component

        self.n_obs += n_obs_minibatch
        self.n_batches += 1

    def compute(self):
        avg_elbo_component = self.elbo_component / self.n_obs

        return avg_elbo_component


class ReconstructionLossMetric(BaseElboMetric):
    """
    Reconstruction loss metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    mode
        Train or validation, used for logging names
    dist_sync_on_step
        optional, by default False
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    def __init__(
        self,
        mode: Literal["train", "validation"],
        dist_sync_on_step: bool = False,
        **kwargs
    ):
        super().__init__(mode=mode, dist_sync_on_step=dist_sync_on_step, **kwargs)

    def update(
        self,
        reconstruction_loss_sum: torch.Tensor,
        n_obs_minibatch: int,
    ):
        """Updates all metrics."""
        super().update(reconstruction_loss_sum, n_obs_minibatch)


class KLLocalMetric(BaseElboMetric):
    """
    KL local loss metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    mode
        Train or validation, used for logging names
    dist_sync_on_step
        optional, by default False
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    def __init__(
        self,
        mode: Literal["train", "validation"],
        dist_sync_on_step: bool = False,
        **kwargs
    ):
        super().__init__(mode=mode, dist_sync_on_step=dist_sync_on_step, **kwargs)

    def update(
        self,
        kl_local_sum: torch.Tensor,
        n_obs_minibatch: int,
    ):
        """Updates all metrics."""
        super().update(kl_local_sum, n_obs_minibatch)


class KLGlobalMetric(BaseElboMetric):
    """
    KL global loss metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    mode
        Train or validation, used for logging names
    dist_sync_on_step
        optional, by default False
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    def __init__(
        self,
        mode: Literal["train", "validation"],
        dist_sync_on_step: bool = False,
        **kwargs
    ):
        super().__init__(mode=mode, dist_sync_on_step=dist_sync_on_step, **kwargs)

    def update(
        self,
        kl_global: torch.Tensor,
        n_obs_minibatch: int,
    ):
        """Updates all metrics."""
        super().update(kl_global, n_obs_minibatch)

    def compute(self):
        avg_elbo_component = self.elbo_component / self.n_batches

        return avg_elbo_component
