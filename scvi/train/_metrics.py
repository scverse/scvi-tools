import torch
from torchmetrics import Metric

from scvi._compat import Literal


class ElboMetric(Metric):
    """
    Elbo metric aggregator for scvi-tools experiments.

    Parameters
    ----------
    n_obs_total
        Number of total observations, for rescaling the ELBO
    mode
        Train or validation, used for logging names
    dist_sync_on_step
        optional, by default False
    **kwargs
        Keyword args for :class:`torchmetrics.Metric`
    """

    def __init__(
        self,
        n_obs_total: int,
        mode: Literal["train", "validation"],
        dist_sync_on_step: bool = False,
        **kwargs
    ):
        super().__init__(dist_sync_on_step=dist_sync_on_step, **kwargs)

        self.n_obs_total = 1 if n_obs_total is None else n_obs_total
        self._mode = mode

        default_val = torch.tensor(0.0)
        self.add_state("reconstruction_loss", default=default_val)
        self.add_state("kl_local", default=default_val)
        self.add_state("kl_global", default=default_val)
        self.add_state("n_obs", default=default_val)
        self.add_state("n_batches", default=default_val)

    @property
    def mode(self):
        return self._mode

    def update(
        self,
        reconstruction_loss_sum: torch.Tensor,
        kl_local_sum: torch.Tensor,
        kl_global: torch.Tensor,
        n_obs_minibatch: int,
    ):
        """Updates all metrics."""
        reconstruction_loss_sum = reconstruction_loss_sum.detach()
        kl_local_sum = kl_local_sum.detach()
        kl_global = kl_global.detach()

        self.reconstruction_loss += reconstruction_loss_sum
        self.kl_local += kl_local_sum
        self.kl_global += kl_global
        self.n_obs += n_obs_minibatch
        self.n_batches += 1

    def compute(self):
        avg_reconstruction_loss = self.reconstruction_loss / self.n_obs
        avg_kl_local = self.kl_local / self.n_obs
        avg_kl_global = self.kl_global / self.n_batches
        # elbo on the scale of one observation
        elbo = (
            avg_reconstruction_loss + avg_kl_local + (avg_kl_global / self.n_obs_total)
        )

        return elbo
