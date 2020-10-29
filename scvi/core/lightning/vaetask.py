from typing import Union

import pytorch_lightning as pl
import torch

from scvi.core.modules._base._base_module import AbstractVAE


class VAETask(pl.LightningModule):
    def __init__(
        self,
        vae_model: AbstractVAE,
        lr=1e-3,
        weight_decay=1e-6,
        n_iter_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
    ):
        super(VAETask, self).__init__()

        self.model = vae_model
        self.lr = lr
        self.weight_decay = weight_decay

    def forward(self, *args, **kwargs):
        """Passthrough to model.forward()."""
        return self.model(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        _, _, scvi_loss = self.forward(batch, kl_weight=self.kl_weight)
        return scvi_loss.loss

    def validation_step(self, batch, batch_idx, optimizer_idx=0):
        _, _, scvi_loss = self.forward(batch)

        loss = scvi_loss.loss
        self.log("elbo_validation", loss)
        return

    def configure_optimizers(self):
        params = filter(lambda p: p.requires_grad, self.model.parameters())
        optimizer = torch.optim.Adam(
            params, lr=self.lr, eps=0.01, weight_decay=self.weight_decay
        )

        return optimizer

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence during training."""
        epoch_criterion = self.n_epochs_kl_warmup is not None
        iter_criterion = self.n_iter_kl_warmup is not None
        if epoch_criterion:
            kl_weight = min(1.0, self.current_epoch / self.n_epochs_kl_warmup)
        elif iter_criterion:
            kl_weight = min(1.0, self.global_step / self.n_iter_kl_warmup)
        else:
            kl_weight = 1.0
        return kl_weight
