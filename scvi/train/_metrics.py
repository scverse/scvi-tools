from typing import List

import torch
from torchmetrics import Metric


class VIMetrics(Metric):
    def __init__(
        self, n_obs_total: int, dist_sync_on_step=False, additional_keys: List = None
    ):
        """Metric aggregator for scvi-tools experiments.

        A number of metrics commonly used in scvi-tools models are tracked down, as well as additional metrics

        Parameters
        ----------
        dist_sync_on_step : bool, optional, by default False
        additional_keys : List, optional
            List of additional metric names to track, by default None
        """
        super().__init__(dist_sync_on_step=dist_sync_on_step)

        if n_obs_total is None:
            self.n_obs_total = 1

        default_val = torch.tensor(0.0)
        self.add_state("reconstruction_loss", default=default_val)
        self.add_state("kl_local", default=default_val)
        self.add_state("kl_global", default=default_val)
        self.add_state("n_obs", default=default_val)
        self.add_state("n_batches", default=default_val)

        self.default_metrics = [
            "elbo",
            "reconstruction_loss",
            "kl_local",
            "elbo",
        ]  # Default keys used for logging
        additional_keys = [] if additional_keys is None else additional_keys
        for new_key in additional_keys:
            self.add_state(new_key, default=default_val)
        self.additional_keys = additional_keys
        self.set_dtype(torch.float)

    def update(
        self,
        reconstruction_loss_sum: torch.Tensor,
        kl_local_sum: torch.Tensor,
        kl_global: torch.Tensor,
        n_obs_minibatch: int,
        **kwargs,
    ):
        """Updates all (or some) metrics."""
        reconstruction_loss_sum = reconstruction_loss_sum.detach().cpu()
        kl_local_sum = kl_local_sum.detach().cpu()
        kl_global = kl_global.cpu()

        self.reconstruction_loss += reconstruction_loss_sum
        self.kl_local += kl_local_sum
        self.kl_global += kl_global
        self.n_obs += n_obs_minibatch
        self.n_batches += 1

        for new_key, tensor_value in kwargs.items():
            old_value = getattr(self, new_key)
            setattr(self, new_key, old_value + tensor_value)

    def compute(self) -> dict:
        avg_reconstruction_loss = self.reconstruction_loss / self.n_obs
        avg_kl_local = self.kl_local / self.n_obs
        kl_global = self.kl_global / self.n_batches
        # elbo on the scale of one observation
        elbo = avg_reconstruction_loss + avg_kl_local + (kl_global / self.n_obs_total)
        main_metrics = {
            "elbo": elbo,
            "reconstruction_loss": avg_reconstruction_loss,
            "kl_local": avg_kl_local,
            "kl_global": kl_global,
        }
        additional_metrics = {key: getattr(self, key) for key in self.additional_keys}
        return {**main_metrics, **additional_metrics}
