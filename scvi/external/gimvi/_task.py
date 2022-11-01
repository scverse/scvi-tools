import torch

from scvi import REGISTRY_KEYS
from scvi.module import Classifier
from scvi.train import METRIC_KEYS, AdversarialTrainingPlan, get_metric_key


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

    def training_step(self, batch, batch_idx, optimizer_idx=0):
        """Training step."""
        kappa = (
            1 - self.kl_weight
            if self.scale_adversarial_loss == "auto"
            else self.scale_adversarial_loss
        )
        if optimizer_idx == 0:
            # batch contains both data loader outputs
            loss_output_objs = []
            n_obs = 0
            zs = []
            for (i, tensors) in enumerate(batch):
                n_obs += tensors[REGISTRY_KEYS.X_KEY].shape[0]
                self.loss_kwargs.update(dict(kl_weight=self.kl_weight, mode=i))
                inference_kwargs = dict(mode=i)
                generative_kwargs = dict(mode=i)
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
                torch.zeros((z.shape[0], 1), device=z.device) + i
                for i, z in enumerate(zs)
            ]
            if kappa > 0 and self.adversarial_classifier is not False:
                fool_loss = self.loss_adversarial_classifier(
                    torch.cat(zs), torch.cat(batch_tensor), False
                )
                loss += fool_loss * kappa

            return {
                METRIC_KEYS.LOSS: loss,
                get_metric_key(
                    METRIC_KEYS.REC_LOSS, metric_type=METRIC_KEYS.SUM
                ): rec_loss,
                get_metric_key(METRIC_KEYS.KL_LOCAL, metric_type=METRIC_KEYS.SUM): kl,
                METRIC_KEYS.KL_GLOBAL: 0.0,
                METRIC_KEYS.N_OBS: n_obs,
            }

        # train adversarial classifier
        # this condition will not be met unless self.adversarial_classifier is not False
        if optimizer_idx == 1:
            zs = []
            for (i, tensors) in enumerate(batch):
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

            return loss

    def validation_step(self, batch, batch_idx, dataloader_idx):
        """Validation step."""
        self.loss_kwargs.update(dict(kl_weight=self.kl_weight, mode=dataloader_idx))
        inference_kwargs = dict(mode=dataloader_idx)
        generative_kwargs = dict(mode=dataloader_idx)
        _, _, loss_output = self.forward(
            batch,
            loss_kwargs=self.loss_kwargs,
            inference_kwargs=inference_kwargs,
            generative_kwargs=generative_kwargs,
        )
        reconstruction_loss = loss_output.reconstruction_loss_sum
        return {
            get_metric_key(
                METRIC_KEYS.REC_LOSS, metric_type=METRIC_KEYS.SUM
            ): reconstruction_loss,
            get_metric_key(
                METRIC_KEYS.KL_LOCAL, metric_type=METRIC_KEYS.SUM
            ): loss_output.kl_local_sum,
            METRIC_KEYS.KL_GLOBAL: loss_output.kl_global,
            METRIC_KEYS.N_OBS: loss_output.n_obs_minibatch,
        }

    def validation_epoch_end(self, outputs):
        """Aggregate validation step information."""
        n_obs, elbo, rec_loss, kl_local = 0, 0, 0, 0
        for dl_out in outputs:
            for tensors in dl_out:
                elbo += (
                    tensors[
                        get_metric_key(
                            METRIC_KEYS.REC_LOSS, metric_type=METRIC_KEYS.SUM
                        )
                    ]
                    + tensors[
                        get_metric_key(
                            METRIC_KEYS.KL_LOCAL, metric_type=METRIC_KEYS.SUM
                        )
                    ]
                )
                rec_loss += tensors[
                    get_metric_key(METRIC_KEYS.REC_LOSS, metric_type=METRIC_KEYS.SUM)
                ]
                kl_local += tensors[
                    get_metric_key(METRIC_KEYS.KL_LOCAL, metric_type=METRIC_KEYS.SUM)
                ]
                n_obs += tensors[METRIC_KEYS.N_OBS]
        # kl global same for each minibatch
        self.log(
            get_metric_key(METRIC_KEYS.ELBO, mode=METRIC_KEYS.VALIDATION), elbo / n_obs
        )
        self.log(
            get_metric_key(METRIC_KEYS.REC_LOSS, mode=METRIC_KEYS.VALIDATION),
            rec_loss / n_obs,
        )
        self.log(
            get_metric_key(METRIC_KEYS.KL_LOCAL, mode=METRIC_KEYS.VALIDATION),
            kl_local / n_obs,
        )
        self.log(
            get_metric_key(METRIC_KEYS.KL_GLOBAL, mode=METRIC_KEYS.VALIDATION), 0.0
        )
