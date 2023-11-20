import torch

from scvi import REGISTRY_KEYS
from scvi.module import Classifier
from scvi.train import AdversarialTrainingPlan


class GIMVITrainingPlan(AdversarialTrainingPlan):
    """gimVI training plan."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if kwargs["adversarial_classifier"] is True:
            self.n_output_classifier = 2
            self.adversarial_classifier = Classifier(
                n_input=self.module.n_latent,
                n_hidden=32,
                n_labels=self.n_output_classifier,
                n_layers=3,
                logits=True,
            )
        else:
            self.adversarial_classifier = kwargs["adversarial_classifier"]
        self.automatic_optimization = False
        self.validation_step_outputs = []

    def training_step(self, batch, batch_idx):
        """Training step."""
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        opts = self.optimizers()
        if not isinstance(opts, list):
            opt1 = opts
            opt2 = None
        else:
            opt1, opt2 = opts
        # batch contains both data loader outputs
        loss_output_objs = []
        n_obs = 0
        zs = []
        for i, tensors in enumerate(batch):
            n_obs += tensors[REGISTRY_KEYS.X_KEY].shape[0]
            self.loss_kwargs.update({"kl_weight": self.kl_weight, "mode": i})
            inference_kwargs = {"mode": i}
            generative_kwargs = {"mode": i}
            inference_outputs, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )
            zs.append(inference_outputs["z"])
            loss_output_objs.append(loss_output)

        loss = sum([scl.loss for scl in loss_output_objs])
        loss /= n_obs
        rec_loss = sum([scl.reconstruction_loss_sum for scl in loss_output_objs])
        kl = sum([scl.kl_local_sum for scl in loss_output_objs])

        # fool classifier if doing adversarial training
        batch_tensor = [
            torch.zeros((z.shape[0], 1), device=z.device) + i for i, z in enumerate(zs)
        ]
        if kappa > 0 and self.adversarial_classifier is not False:
            fool_loss = self.loss_adversarial_classifier(
                torch.cat(zs), torch.cat(batch_tensor), False
            )
            loss += fool_loss * kappa
        opt1.zero_grad()
        self.manual_backward(loss)
        opt1.step()
        return_dict = {
            "loss": loss,
            "reconstruction_loss_sum": rec_loss,
            "kl_local_sum": kl,
            "kl_global": 0.0,
            "n_obs": n_obs,
        }

        # train adversarial classifier
        # this condition will not be met unless self.adversarial_classifier is not False
        if opt2 is not None:
            zs = []
            for i, tensors in enumerate(batch):
                inference_inputs = self.module._get_inference_input(tensors)
                inference_inputs.update({"mode": i})
                outputs = self.module.inference(**inference_inputs)
                zs.append(outputs["z"])

            batch_tensor = [
                torch.zeros((z.shape[0], 1), device=z.device) + i
                for i, z in enumerate(zs)
            ]
            loss = self.loss_adversarial_classifier(
                torch.cat(zs).detach(), torch.cat(batch_tensor), True
            )
            loss *= kappa
            opt2.zero_grad()
            self.manual_backward(loss)
            opt2.step()

        return return_dict

    def validation_step(self, batch, batch_idx, dataloader_idx):
        """Validation step."""
        self.loss_kwargs.update({"kl_weight": self.kl_weight, "mode": dataloader_idx})
        inference_kwargs = {"mode": dataloader_idx}
        generative_kwargs = {"mode": dataloader_idx}
        _, _, loss_output = self.forward(
            batch,
            loss_kwargs=self.loss_kwargs,
            inference_kwargs=inference_kwargs,
            generative_kwargs=generative_kwargs,
        )
        reconstruction_loss = loss_output.reconstruction_loss_sum
        self.validation_step_outputs.append(
            {
                "reconstruction_loss_sum": reconstruction_loss,
                "kl_local_sum": loss_output.kl_local_sum,
                "kl_global": loss_output.kl_global,
                "n_obs": loss_output.n_obs_minibatch,
            }
        )

    def on_validation_epoch_end(self):
        """Aggregate validation step information."""
        super().on_validation_epoch_end()
        outputs = self.validation_step_outputs
        n_obs, elbo, rec_loss, kl_local = 0, 0, 0, 0
        for val_metrics in outputs:
            elbo += val_metrics["reconstruction_loss_sum"] + val_metrics["kl_local_sum"]
            rec_loss += val_metrics["reconstruction_loss_sum"]
            kl_local += val_metrics["kl_local_sum"]
            n_obs += val_metrics["n_obs"]
        # kl global same for each minibatch
        self.log("elbo_validation", elbo / n_obs)
        self.log("reconstruction_loss_validation", rec_loss / n_obs)
        self.log("kl_local_validation", kl_local / n_obs)
        self.log("kl_global_validation", 0.0)
        self.validation_step_outputs.clear()  # free memory
