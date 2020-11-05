from typing import Union

import pytorch_lightning as pl
import torch
from torch.nn import functional as F

from scvi.core.modules._base._base_module import AbstractVAE
from scvi import _CONSTANTS


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
        self.n_iter_kl_warmup = n_iter_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup

    def forward(self, *args, **kwargs):
        """Passthrough to model.forward()."""
        return self.model(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        loss_kwargs = dict(kl_weight=self.kl_weight)
        _, _, scvi_loss = self.forward(batch, loss_kwargs=loss_kwargs)
        return scvi_loss.loss

    # on validation epoch end
    # mean of the mean etc
    # should return scvi_loss object
    # on_validation_epoch_end should iterate over it
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


class ClassifierTask(VAETask):
    def __init__(self, vae_model: AbstractVAE):
        self.model = vae_model

    def training_step(self, batch, batch_idx):
        x = batch[_CONSTANTS.X_KEY]
        b = batch[_CONSTANTS.BATCH_KEY]
        labels_train = batch[_CONSTANTS.LABELS_KEY]
        if self.sampling_model:
            # TODO: we need to document that it looks for classify
            if hasattr(self.sampling_model, "classify"):
                return F.cross_entropy(
                    self.sampling_model.classify(x, b), labels_train.view(-1)
                )
            else:
                if self.sampling_model.log_variational:
                    x = torch.log(1 + x)
                if self.sampling_zl:
                    x_z = self.sampling_model.z_encoder(x, b)[0]
                    x_l = self.sampling_model.l_encoder(x, b)[0]
                    x = torch.cat((x_z, x_l), dim=-1)
                else:
                    x = self.sampling_model.z_encoder(x, b)[0]
        return F.cross_entropy(self.model(x), labels_train.view(-1))


class SemiSupervisedTask(VAETask):
    def __init__(
        self,
        vae_model: AbstractVAE,
        lr=1e-3,
        weight_decay=1e-6,
        n_iter_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
    ):
        super(SemiSupervisedTask, self).__init__()
        self.classifier_trainer = ClassifierTask()

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        full_dataset = batch[0]
        labelled_dataset = batch[1]
        input_kwargs = dict(feed_labels=False)
        _, _, scvi_losses = self.forward(full_dataset, loss_kwargs=input_kwargs)
        loss = scvi_losses.loss
        x = labelled_dataset[_CONSTANTS.X_KEY]
        y = labelled_dataset[_CONSTANTS.LABELS_KEY]
        classification_loss = F.cross_entropy(self.model.classify(x), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss

    def on_epoch_end(self):
        self.model.eval()
        self.classifier_trainer.train(
            self.n_epochs_classifier, lr=self.lr_classification
        )
        self.model.train()
        return super().on_epoch_end()
