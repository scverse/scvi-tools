from typing import Union

import pytorch_lightning as pl
import torch
from torch.nn import functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau

from scvi._compat import Literal
from scvi.core.modules import Classifier
from scvi.core.modules._base._base_module import AbstractVAE
from scvi.core.modules.utils import one_hot
from scvi import _CONSTANTS


class VAETask(pl.LightningModule):
    """
    Lightning module task to train scvi-tools modules

    Parameters
    ----------
    vae_model
        A model instance from class ``AbstractVAE``
    """

    def __init__(
        self,
        vae_model: AbstractVAE,
        lr=1e-3,
        weight_decay=1e-6,
        n_iter_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        adversarial_classifier: Union[bool, Classifier] = False,
        scale_adversarial_loss: Union[float, Literal["auto"]] = "auto",
    ):
        super(VAETask, self).__init__()
        self.model = vae_model
        self.lr = lr
        self.weight_decay = weight_decay
        self.n_iter_kl_warmup = n_iter_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.reduce_lr_on_plateau = reduce_lr_on_plateau
        self.lr_factor = lr_factor
        self.lr_patience = lr_patience
        self.lr_scheduler_metric = lr_scheduler_metric

        if adversarial_classifier is True:
            self.adversarial_classifier = Classifier(
                n_input=self.model.n_latent,
                n_hidden=32,
                n_labels=self.model.n_batch,
                n_layers=2,
                logits=True,
            )
        else:
            self.adversarial_classifier = adversarial_classifier
        self.scale_adversarial_loss = scale_adversarial_loss

    def forward(self, *args, **kwargs):
        """Passthrough to model.forward()."""
        return self.model(*args, **kwargs)

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        batch_tensor = batch[_CONSTANTS.BATCH_KEY]
        if optimizer_idx == 0:
            loss_kwargs = dict(kl_weight=self.kl_weight)
            inference_outputs, _, scvi_loss = self.forward(
                batch, loss_kwargs=loss_kwargs
            )
            # fool classifier if doing adversarial training
            if kappa > 0 and self.adversarial_classifier is not False:
                z = inference_outputs["z"]
                fool_loss = self.loss_discriminator(z, batch_tensor, False)
                scvi_loss.loss += fool_loss * kappa
            return scvi_loss.loss

        # train adversarial classifier
        # this condition will not be met unless self.adversarial_classifier is not False
        if optimizer_idx == 1:
            inference_inputs = self.model._get_inference_input(batch)
            outputs = self.model.inference(**inference_inputs)
            z = outputs["z"]
            adversarial_loss = self.loss_discriminator(z.detach(), batch_tensor, True)
            adversarial_loss *= kappa
            return adversarial_loss

    def validation_step(self, batch, batch_idx):
        _, _, scvi_loss = self.forward(batch)
        reconstruction_loss = scvi_loss.reconstruction_loss
        return {
            "reconstruction_loss_sum": reconstruction_loss.sum(),
            "kl_local_sum": scvi_loss.kl_local.sum(),
            "kl_global": scvi_loss.kl_global,
            "n_obs": reconstruction_loss.shape[0],
        }

    def validation_epoch_end(self, outputs):
        """Aggregate validation step information."""
        n_obs, elbo, rec_loss, kl_local = 0, 0, 0, 0
        for tensors in outputs:
            elbo += tensors["reconstruction_loss_sum"] + tensors["kl_local_sum"]
            rec_loss += tensors["reconstruction_loss_sum"]
            kl_local += tensors["kl_local_sum"]
            n_obs += tensors["n_obs"]
        # kl global same for each minibatch
        kl_global = tensors["kl_global"]
        elbo += kl_global

        self.log("elbo_validation", elbo / n_obs)
        self.log("reconstruction_loss_validation", rec_loss / n_obs)
        self.log("kl_local_validation", kl_local / n_obs)
        self.log("kl_global_validation", kl_global)

    def configure_optimizers(self):
        params1 = filter(lambda p: p.requires_grad, self.model.parameters())
        optimizer1 = torch.optim.Adam(
            params1, lr=self.lr, eps=0.01, weight_decay=self.weight_decay
        )
        config1 = {"optimizer": optimizer1}
        if self.reduce_lr_on_plateau:
            scheduler1 = ReduceLROnPlateau(
                optimizer1, patience=self.lr_patience, factor=self.lr_factor
            )
            config1.update(
                {
                    "lr_scheduler": scheduler1,
                    "monitor": self.lr_scheduler_metric,
                },
            )

        if self.adversarial_classifier is not False:
            params2 = filter(
                lambda p: p.requires_grad, self.adversarial_classifier.parameters()
            )
            optimizer2 = torch.optim.Adam(
                params2, lr=1e-3, eps=0.01, weight_decay=self.weight_decay
            )
            config2 = {"optimizer": optimizer2}
            return config1, config2

        return config1

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

    def loss_discriminator(self, z, batch_index, predict_true_class=True):
        n_classes = self.model.n_batch
        cls_logits = torch.nn.LogSoftmax(dim=1)(self.discriminator(z))

        if predict_true_class:
            cls_target = one_hot(batch_index, n_classes)
        else:
            one_hot_batch = one_hot(batch_index, n_classes)
            cls_target = torch.zeros_like(one_hot_batch)
            # place zeroes where true label is
            cls_target.masked_scatter_(
                ~one_hot_batch.bool(), torch.ones_like(one_hot_batch) / (n_classes - 1)
            )

        l_soft = cls_logits * cls_target
        loss = -l_soft.sum(dim=1).mean()

        return loss


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
