import torch
from torchmetrics import Metric


class ElboMetric(Metric):
    def __init__(self, n_obs_total: int, dist_sync_on_step: bool = False):
        """
        Elbo metric aggregator for scvi-tools experiments.

        Parameters
        ----------
        n_obs_total
            Number of total observations, for rescaling the ELBO
        dist_sync_on_step
            optional, by default False
        """
        super().__init__(dist_sync_on_step=dist_sync_on_step)

        self.n_obs_total = 1 if n_obs_total is None else n_obs_total

        default_val = torch.tensor(0.0)
        self.add_state("reconstruction_loss", default=default_val)
        self.add_state("kl_local", default=default_val)
        self.add_state("kl_global", default=default_val)
        self.add_state("n_obs", default=default_val)
        self.add_state("n_batches", default=default_val)

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
        kl_global = kl_global

        self.reconstruction_loss += reconstruction_loss_sum
        self.kl_local += kl_local_sum
        self.kl_global += kl_global
        self.n_obs += n_obs_minibatch
        self.n_batches += 1

    def compute(self):
        avg_reconstruction_loss = self.reconstruction_loss / self.n_obs
        avg_kl_local = self.kl_local / self.n_obs
        kl_global = self.kl_global / self.n_batches
        # elbo on the scale of one observation
        elbo = avg_reconstruction_loss + avg_kl_local + (kl_global / self.n_obs_total)

        return elbo
