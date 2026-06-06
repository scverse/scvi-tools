from __future__ import annotations

import torch

from scvi.train import TrainingPlan


class SCPOLITrainingPlan(TrainingPlan):
    """
    Training plan for scPoli with two-phase training.

    Parameters
    ----------
    module
        The SCPOLIModule instance.
    pretraining_epochs
        Number of epochs to train with prototype loss only (no ELBO).
    eta
        Weight for prototype loss during full training.
    lr
        Learning rate.
    weight_decay
        Weight decay for optimizer.
    """

    def __init__(
        self,
        module,
        pretraining_epochs: int = 50,
        eta: float = 5.0,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        **kwargs,
    ):
        super().__init__(module, lr=lr, weight_decay=weight_decay, **kwargs)
        self.pretraining_epochs = pretraining_epochs
        self.eta = eta

    def training_step(self, batch, batch_idx):
        """Two-phase training step."""
        # --- Phase 1: Pretraining ---
        # Only prototype loss is active
        # Encoder learns to cluster cells by cell type
        if self.current_epoch < self.pretraining_epochs:
            loss_output = self.module.loss(
                tensors=batch,
                inference_outputs=self.module.inference(**self.module._get_inference_input(batch)),
                generative_outputs=self.module.generative(
                    **self.module._get_generative_input(
                        batch, self.module.inference(**self.module._get_inference_input(batch))
                    )
                ),
                kl_weight=0.0,  # no KL during pretraining
                proto_weight=self.eta,
            )

        # --- Phase 2: Full Training ---
        # ELBO + prototype loss together
        else:
            # KL annealing: gradually increase kl_weight from 0 to 1
            kl_weight = min(1.0, (self.current_epoch - self.pretraining_epochs) / 50)
            loss_output = self.module.loss(
                tensors=batch,
                inference_outputs=self.module.inference(**self.module._get_inference_input(batch)),
                generative_outputs=self.module.generative(
                    **self.module._get_generative_input(
                        batch, self.module.inference(**self.module._get_inference_input(batch))
                    )
                ),
                kl_weight=kl_weight,
                proto_weight=self.eta,
            )

        return loss_output.loss

    def validation_step(self, batch, batch_idx):
        """Validation step — always uses full loss."""
        loss_output = self.module.loss(
            tensors=batch,
            inference_outputs=self.module.inference(**self.module._get_inference_input(batch)),
            generative_outputs=self.module.generative(
                **self.module._get_generative_input(
                    batch, self.module.inference(**self.module._get_inference_input(batch))
                )
            ),
            kl_weight=1.0,
            proto_weight=self.eta,
        )
        self.log("validation_loss", loss_output.loss, on_epoch=True)
        return loss_output.loss

    def configure_optimizers(self):
        """Single Adam optimizer for all parameters."""
        optimizer = torch.optim.Adam(
            self.module.parameters(),
            lr=self.lr,
            weight_decay=self.weight_decay,
        )
        return optimizer
