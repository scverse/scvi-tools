from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau

from scvi import settings
from scvi.module import Classifier
from scvi.train import AdversarialTrainingPlan

if TYPE_CHECKING:
    from typing import Literal

    from scvi.module.base import BaseModuleClass
    from scvi.train._trainingplans import TorchOptimizerCreator

ASSAY_KEY = "assay"
ADVERSARIAL_GROUP_KEY = "adversarial_group"


class SCVIXTrainingPlan(AdversarialTrainingPlan):
    """Adversarial training plan with scVI-X-specific classifier inputs."""

    def __init__(
        self,
        module: BaseModuleClass,
        *,
        optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam",
        optimizer_creator: TorchOptimizerCreator | None = None,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        n_steps_kl_warmup: int = None,
        n_epochs_kl_warmup: int = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        lr_min: float = 0,
        adversarial_classifier: bool | Classifier = False,
        adversarial_key: str = ASSAY_KEY,
        adversarial_steps: int = 1,
        scale_adversarial_loss: float | Literal["auto"] = "auto",
        compile: bool = False,
        compile_kwargs: dict | None = None,
        **loss_kwargs,
    ):
        super().__init__(
            module=module,
            optimizer=optimizer,
            optimizer_creator=optimizer_creator,
            lr=lr,
            weight_decay=weight_decay,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_threshold=lr_threshold,
            lr_scheduler_metric=lr_scheduler_metric,
            lr_min=lr_min,
            adversarial_classifier=False,
            scale_adversarial_loss=scale_adversarial_loss,
            compile=compile,
            compile_kwargs=compile_kwargs,
            **loss_kwargs,
        )
        self.adversarial_key = adversarial_key
        self.adversarial_steps = adversarial_steps

        if adversarial_classifier is True:
            n_adversarial_name = f"n_{adversarial_key}"
            if not hasattr(self.module, n_adversarial_name):
                raise ValueError(
                    f"Adversarial key {adversarial_key!r} not found in module setup args."
                )
            self.n_output_classifier = getattr(self.module, n_adversarial_name)
            if self.n_output_classifier == 1:
                warnings.warn(
                    "Disabling adversarial classifier as there is only one class.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
                self.adversarial_classifier = False
            else:
                self.adversarial_classifier = Classifier(
                    n_input=self.module.n_latent + getattr(self.module, "n_adversarial_group", 0),
                    n_hidden=128,
                    n_labels=self.n_output_classifier,
                    n_layers=2,
                    logits=True,
                    use_batch_norm=False,
                    use_layer_norm=True,
                )
        else:
            self.adversarial_classifier = adversarial_classifier

    def loss_adversarial_classifier(
        self,
        z: torch.Tensor,
        adversarial_group: torch.Tensor,
        class_index: torch.Tensor,
        predict_true_class: bool = True,
    ) -> torch.Tensor:
        """Loss for the scVI-X adversarial classifier."""
        n_classes = self.n_output_classifier
        n_adv_group = getattr(self.module, "n_adversarial_group", 0)
        if n_adv_group > 0:
            adversarial_group_one_hot = torch.nn.functional.one_hot(
                adversarial_group.squeeze(-1).long(), num_classes=n_adv_group
            ).float()
            z = torch.cat([z, adversarial_group_one_hot], dim=1)
        if predict_true_class:
            z = z.detach()

        cls_logits = self.adversarial_classifier(z)
        if predict_true_class:
            return torch.nn.functional.cross_entropy(cls_logits, class_index.squeeze(-1).long())

        one_hot_class = torch.nn.functional.one_hot(class_index.squeeze(-1), n_classes).float()
        cls_target = (1 - one_hot_class) / (n_classes - 1)
        return -(cls_target * torch.nn.functional.log_softmax(cls_logits, dim=1)).sum(dim=1).mean()

    def training_step(self, batch, batch_idx):
        """Training step for scVI-X adversarial training."""
        if "kl_weight" in self.loss_kwargs:
            self.loss_kwargs.update({"kl_weight": self.kl_weight})
            self.log("kl_weight", self.kl_weight, on_step=True, on_epoch=False)
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        class_tensor = batch[self.adversarial_key].long()

        opts = self.optimizers()
        if not isinstance(opts, list):
            opt1 = opts
            opt2 = None
        else:
            opt1, opt2 = opts

        inference_outputs, _, scvi_loss = self.forward(batch, loss_kwargs=self.loss_kwargs)
        z = inference_outputs["z"]
        adversarial_group = inference_outputs.get(ADVERSARIAL_GROUP_KEY)
        if adversarial_group is None:
            adversarial_group = torch.zeros(z.size(0), device=z.device, dtype=torch.long)

        loss = scvi_loss.loss
        orig_loss = loss
        if kappa > 0 and self.adversarial_classifier is not False:
            fool_loss = self.loss_adversarial_classifier(
                z,
                adversarial_group,
                class_tensor,
                predict_true_class=False,
            )
            loss += fool_loss * kappa

        self.log("train_loss", loss, on_step=self.on_step, on_epoch=self.on_epoch, prog_bar=True)
        if self.on_step:
            self.trainer.logger.log_metrics(
                {"train_loss_step": loss},
                step=self.global_step,
            )
        self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")
        opt1.zero_grad()
        self.manual_backward(loss)
        opt1.step()

        if opt2 is not None and kappa > 0:
            qz = inference_outputs["qz"]
            for _ in range(max(self.adversarial_steps, 0)):
                z_classifier = qz.sample()
                classifier_loss = kappa * self.loss_adversarial_classifier(
                    z_classifier,
                    adversarial_group,
                    class_tensor,
                    predict_true_class=True,
                )
                opt2.zero_grad()
                self.manual_backward(classifier_loss)
                opt2.step()

        if scvi_loss.extra_metrics is not None and len(scvi_loss.extra_metrics.keys()) > 0:
            self.prepare_scib_autotune(scvi_loss.extra_metrics, "training")

        return orig_loss

    def configure_optimizers(self):
        """Configure optimizers for scVI-X adversarial training."""
        params1 = filter(lambda p: p.requires_grad, self.module.parameters())
        optimizer1 = self.get_optimizer_creator()(params1)
        config1 = {"optimizer": optimizer1}
        if self.reduce_lr_on_plateau:
            scheduler1 = ReduceLROnPlateau(
                optimizer1,
                patience=self.lr_patience,
                factor=self.lr_factor,
                threshold=self.lr_threshold,
                min_lr=self.lr_min,
                threshold_mode="abs",
            )
            config1.update(
                {
                    "lr_scheduler": {
                        "scheduler": scheduler1,
                        "monitor": self.lr_scheduler_metric,
                    },
                },
            )

        if self.adversarial_classifier is not False:
            params2 = filter(lambda p: p.requires_grad, self.adversarial_classifier.parameters())
            optimizer2 = torch.optim.Adam(params2, lr=3e-4, eps=1e-4, weight_decay=1e-9)
            opts = [config1.pop("optimizer"), optimizer2]
            if "lr_scheduler" in config1:
                return opts, [config1["lr_scheduler"]]
            return opts

        return config1
