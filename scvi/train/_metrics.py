from typing import List

import torch
from torchmetrics import Metric


class VIMetrics(Metric):
    def __init__(self, dist_sync_on_step=False, additional_keys: List = None):
        """Metric aggregator for scvi-tools experiments.

        A number of metrics commonly used in scvi-tools models are tracked down, as well as additional metrics

        Parameters
        ----------
        dist_sync_on_step : bool, optional, by default False
        additional_keys : List, optional
            List of additional metric names to track, by default None
        """
        super().__init__(dist_sync_on_step=dist_sync_on_step)

        default_val = torch.tensor(0)
        # self.add_state("elbo", default=default_val)
        self.add_state("reconstruction_loss", default=default_val)
        self.add_state("kl_local", default=default_val)
        self.add_state("kl_global", default=default_val)
        self.add_state("n_obs", default=default_val)

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
        scvi_losses,
        **kwargs,
    ):
        """Updates all (or some) metrics."""
        n_obs = scvi_losses.reconstruction_loss.shape[0]
        kl_local_sum = scvi_losses.kl_local.sum().detach().cpu()
        kl_global = scvi_losses.kl_global.cpu()
        reconstruction_loss_sum = scvi_losses.reconstruction_loss.sum().detach().cpu()

        self.reconstruction_loss += reconstruction_loss_sum
        self.kl_local += kl_local_sum
        self.kl_global += kl_global
        self.n_obs += n_obs

        for new_key, tensor_value in kwargs.items():
            old_value = getattr(self, new_key)
            setattr(self, new_key, old_value + tensor_value)

    def compute(self):
        reconstruction_loss = self.reconstruction_loss.squeeze() / self.n_obs
        kl_local = self.kl_local.squeeze() / self.n_obs
        kl_global = (
            self.kl_global / self.n_obs
        ).squeeze()  # Taking mean of kl global accross batches
        elbo = reconstruction_loss + kl_local + (kl_global / self.n_obs)
        main_metrics = {
            "elbo": elbo,
            "reconstruction_loss": reconstruction_loss,
            "kl_local": kl_local,
            "kl_global": kl_global,
        }
        additional_metrics = {
            key: getattr(self, key).squeeze() for key in self.additional_keys
        }
        return {**main_metrics, **additional_metrics}
