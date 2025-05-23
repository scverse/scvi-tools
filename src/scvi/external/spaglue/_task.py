import torch

from scvi import REGISTRY_KEYS
from scvi.train import TrainingPlan


class SPAGLUETrainingPlan(TrainingPlan):
    def __init__(self, *args: tuple, **kwargs: dict) -> None:
        super().__init__(*args, **kwargs)

        self.history = {
            "train_loss": [],  # Batch-level training losses
            "validation_loss": [],  # Batch-level validation losses
            "train_loss_rec": [],  # Batch-level training losses
            "validation_loss_rec": [],  # Batch-level validation losses
            "train_loss_kl": [],  # Batch-level training losses
            "validation_loss_kl": [],  # Batch-level validation losses
            "avg_train_loss": [],  # Epoch-level average training losses
            "avg_validation_loss": [],  # Epoch-level average validation losses
            "avg_train_loss_rec": [],  # Epoch-level average training losses
            "avg_validation_loss_rec": [],  # Epoch-level average validation losses
            "avg_train_loss_kl": [],  # Epoch-level average training losses
            "avg_validation_loss_kl": [],  # Epoch-level average validation losses
            "avg_train_loss_seq": [],
            "train_loss_seq": [],
            "avg_train_loss_spa": [],
            "train_loss_spa": [],
            "bias_params_diss": [],
            "dispersion_params_diss": [],
            "bias_params_spa": [],
            "dispersion_params_spa": [],
        }

        self.automatic_optimization = False  # important for adversarial setup

    def training_step(self, batch: list[dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
        # batch contains output of both dataloaders
        # There is a possibility to give batch_idx argument
        # (can be used e.g. for gradient accumulation)
        """Training step."""
        opt = self.optimizers()

        # batch contains both data loader outputs (list of 2 train_dl)
        loss_output_objs = []
        # n_obs = 0
        for i, tensors in enumerate(
            batch
        ):  # contains data from all datasets - mode is either 0 or 1
            n_obs = tensors[REGISTRY_KEYS.X_KEY].shape[0]
            n_var = tensors[REGISTRY_KEYS.X_KEY].shape[1]
            self.loss_kwargs.update({"mode": i})
            inference_kwargs = {"mode": i}
            generative_kwargs = {"mode": i}
            _, _, loss_output = self.forward(  # perform inf and gen steps
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)

            kl_divergence = loss_output.kl_local["kl_local"]
            kl_divergence = torch.sum(kl_divergence) / (n_obs * n_var)

            loss = loss_output.loss

            loss_dict = {
                "reconstruction_loss": reconstruction_loss,
                "kl_divergence": kl_divergence,
                "loss": loss,
            }

            loss_output_objs.append(loss_dict)

        loss = sum(i["loss"] for i in loss_output_objs) / 2
        rec_loss = sum(i["reconstruction_loss"] for i in loss_output_objs) / 2
        kl_loss = sum(i["kl_divergence"] for i in loss_output_objs) / 2

        self.log("train_loss", loss, prog_bar=True, on_step=True)

        self.history["train_loss"].append(loss.item())
        self.history["train_loss_rec"].append(rec_loss.item())
        self.history["train_loss_kl"].append(kl_loss.item())

        opt.zero_grad()  # zero out gradient
        self.manual_backward(loss)  # recompute gradients via backpropagation
        opt.step()  # updates the model weights

        return {"loss": loss}

    def validation_step(self, batch: list[dict[str, torch.Tensor]]) -> None:
        """Validation step."""
        loss_output_objs = []

        for i, tensors in enumerate(
            batch
        ):  # contains data from all datasets - mode is either 0 or 1
            n_obs = tensors[REGISTRY_KEYS.X_KEY].shape[0]
            n_var = tensors[REGISTRY_KEYS.X_KEY].shape[1]

            self.loss_kwargs.update({"mode": i})
            inference_kwargs = {"mode": i}
            generative_kwargs = {"mode": i}
            _, _, loss_output = self.forward(  # perform inf and gen steps
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )
            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)

            kl_divergence = loss_output.kl_local["kl_local"]
            kl_divergence = torch.sum(kl_divergence) / (n_obs * n_var)

            loss = loss_output.loss

            loss_dict = {
                "reconstruction_loss": reconstruction_loss,
                "kl_divergence": kl_divergence,
                "loss": loss,
            }

            loss_output_objs.append(loss_dict)

            # look at the kl divergences for both modalities separately
            if i == 0:
                self.history["train_loss_seq"].append(kl_divergence.item())
            if i == 1:
                self.history["train_loss_spa"].append(kl_divergence.item())

        loss = sum(i["loss"] for i in loss_output_objs) / 2
        rec_loss = sum(i["reconstruction_loss"] for i in loss_output_objs) / 2
        kl_loss = sum(i["kl_divergence"] for i in loss_output_objs) / 2

        self.log("validation_loss", loss, prog_bar=True, on_step=True)

        self.history["validation_loss"].append(loss.item())
        self.history["validation_loss_rec"].append(rec_loss.item())
        self.history["validation_loss_kl"].append(kl_loss.item())

    def on_train_epoch_end(self) -> None:
        """Log training loss at the end of each epoch."""
        if len(self.history["train_loss"]) > 0:
            avg_train_loss = sum(self.history["train_loss"]) / len(self.history["train_loss"])
            avg_train_loss_rec = sum(self.history["train_loss_rec"]) / len(
                self.history["train_loss_rec"]
            )
            avg_train_loss_kl = sum(self.history["train_loss_kl"]) / len(
                self.history["train_loss_kl"]
            )
            avg_train_loss_seq = sum(self.history["train_loss_seq"]) / len(
                self.history["train_loss_seq"]
            )
            avg_train_loss_spa = sum(self.history["train_loss_spa"]) / len(
                self.history["train_loss_spa"]
            )

            self.history["avg_train_loss"].append(avg_train_loss)
            self.history["avg_train_loss_rec"].append(avg_train_loss_rec)
            self.history["avg_train_loss_kl"].append(avg_train_loss_kl)
            self.history["avg_train_loss_seq"].append(avg_train_loss_seq)
            self.history["avg_train_loss_spa"].append(avg_train_loss_spa)

            self.log("avg_train_loss", avg_train_loss, prog_bar=True)

        # Extract and log model parameters
        bias_params_diss = self.module.z_decoder_diss.bias.detach().cpu().numpy()
        dispersion_params_diss = self.module.z_decoder_diss.log_theta.exp().detach().cpu().numpy()

        bias_params_spa = self.module.z_decoder_spa.bias.detach().cpu().numpy()
        dispersion_params_spa = self.module.z_decoder_spa.log_theta.exp().detach().cpu().numpy()

        self.history["bias_params_diss"].append(bias_params_diss)
        self.history["dispersion_params_diss"].append(dispersion_params_diss)
        self.history["bias_params_spa"].append(bias_params_spa)
        self.history["dispersion_params_spa"].append(dispersion_params_spa)

        self.history["train_loss"].clear()
        self.history["train_loss_rec"].clear()
        self.history["train_loss_kl"].clear()
        self.history["train_loss_seq"].clear()
        self.history["train_loss_spa"].clear()

    def on_validation_epoch_end(self) -> None:
        """Log validation loss at the end of each epoch."""
        if len(self.history["validation_loss"]) > 0:
            avg_val_loss = sum(self.history["validation_loss"]) / len(
                self.history["validation_loss"]
            )
            avg_val_loss_rec = sum(self.history["validation_loss_rec"]) / len(
                self.history["validation_loss_rec"]
            )
            avg_val_loss_kl = sum(self.history["validation_loss_kl"]) / len(
                self.history["validation_loss_kl"]
            )

            self.history["avg_validation_loss"].append(avg_val_loss)
            self.history["avg_validation_loss_rec"].append(avg_val_loss_rec)
            self.history["avg_validation_loss_kl"].append(avg_val_loss_kl)

            # Log epoch-level losses
            self.log("avg_validation_loss", avg_val_loss, prog_bar=True)
            self.log("avg_validation_loss_rec", avg_val_loss_rec, prog_bar=False)
            self.log("avg_validation_loss_kl", avg_val_loss_kl, prog_bar=False)

        self.history["validation_loss"].clear()
        self.history["validation_loss_rec"].clear()
        self.history["validation_loss_kl"].clear()
