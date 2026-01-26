import torch

from scvi import REGISTRY_KEYS
from scvi.external.diagvi._utils import (
    compute_graph_loss,
    kl_divergence_graph,
)
from scvi.train import TrainingPlan
from scvi.utils import dependencies


def _anneal_param(
    current_epoch: int,
    max_epochs: int,
    init_value: float,
    target_value: float,
) -> float:
    """Linearly anneal a parameter from init_value to target_value over one third of max_epochs.

    Parameters
    ----------
    current_epoch
        Current epoch number.
    max_epochs
        Total number of epochs.
    init_value
        Initial value of the parameter.
    target_value
        Target value of the parameter.

    Returns
    -------
    float
        Annealed parameter value.
    """
    anneal_epochs = max_epochs // 3
    if current_epoch >= anneal_epochs:
        return target_value
    else:
        progress = current_epoch / anneal_epochs
        return init_value + (target_value - init_value) * progress


class DiagTrainingPlan(TrainingPlan):
    """Training plan for DiagVI model.

    Parameters
    ----------
    module
        The DiagVI module to be trained.
    lam_graph
        Weight for the graph loss component.
    lam_kl
        Weight for the KL divergence loss component.
    lam_data
        Weight for the data reconstruction loss component.
    lam_sinkhorn
        Weight for the Sinkhorn loss component.
    lam_class
        Weight for the classification loss component. 
        If None, it will be set dynamically based on the presence of labels.
    sinkhorn_p
        The p parameter for the Sinkhorn loss.
    sinkhorn_blur
        The blur parameter for the Sinkhorn loss.
    sinkhorn_reach
        The reach parameter for the Sinkhorn loss.
    lr
        Learning rate for the optimizer.
    n_epochs_sinkhorn_warmup
        Number of epochs for warming up the Sinkhorn loss weight.
    loss_annealing
        Whether to anneal the Sinkhorn loss parameters over training.
    *args
        Additional positional arguments passed to :class:`~scvi.train.TrainingPlan`.
    **kwargs
        Additional keyword arguments passed to :class:`~scvi.train.TrainingPlan`.
    """
    def __init__(
        self,
        module: torch.nn.Module,
        lam_graph: float = 1.0,
        lam_kl: float = 1.0,
        lam_data: float = 1.0,
        lam_sinkhorn: float = 1.0,
        lam_class: float | None = None,
        sinkhorn_p: int = 2,
        sinkhorn_blur: float = 1.0,
        sinkhorn_reach: float = 1.0,
        lr: float = 1e-3,
        n_epochs_sinkhorn_warmup: int | None = None,
        loss_annealing: bool = True,
        *args,
        **kwargs,
    ) -> None:
        super().__init__(module, *args, **kwargs)

        self.lam_graph = lam_graph
        self.lam_kl = lam_kl
        self.lam_data = lam_data
        self.lam_sinkhorn = lam_sinkhorn
        self.lam_class = lam_class
        self.sinkhorn_p = sinkhorn_p
        self.sinkhorn_reach = sinkhorn_reach
        self.sinkhorn_blur = sinkhorn_blur
        self.lr = lr
        self.n_epochs_sinkhorn_warmup = n_epochs_sinkhorn_warmup
        self.loss_annealing = loss_annealing
        
        # Use larger values initially to do annealing for OT parameters
        self.init_blur = 10 * self.sinkhorn_blur
        self.init_reach = 10 * self.sinkhorn_reach

    @dependencies("geomloss")
    def training_step(self, batch: dict[str, dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
        """Training step.
        
        During training, computes the losses for each modality (NLL, KL, and classification), 
        the graph loss, the Sinkhorn loss between modalities, and combines them into a total loss.

        Parameters
        ----------
        batch
            A batch of data containing tensors for each modality.

        Returns
        -------
        dict[str, torch.Tensor]
            A dictionary containing the total loss for backpropagation.
        """
        import geomloss

        # Compute losses for each modality
        loss_output_objs = []
        for _i, (name, tensors) in enumerate(batch.items()):
            batch_size = tensors[REGISTRY_KEYS.X_KEY].shape[0]

            # Update loss kwargs for current modality
            self.loss_kwargs.update(
                {"lam_kl": self.lam_kl, "lam_data": self.lam_data, "mode": name}
            )
            inference_kwargs = {"mode": name}
            generative_kwargs = {"mode": name}

            # Calculate reconstruction, KL, and classification losses
            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            # Log reconstruction loss
            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)
            self.log(
                f"nll_{name}",
                reconstruction_loss,
                batch_size=batch_size,
                on_epoch=True,
                on_step=False,
            )

            # Log KL divergence loss
            kl_divergence = loss_output.kl_local["kl_local"]
            self.log(
                f"kl_{name}",
                kl_divergence,
                batch_size=batch_size,
                on_epoch=True,
                on_step=False,
            )

            # Log total loss
            loss = loss_output.loss
            self.log(
                f"train_loss_{name}",
                loss,
                batch_size=batch_size,
                on_epoch=True,
                on_step=True,
            )

            loss_dict = {
                "z": loss_output.extra_metrics["z"],
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
                "classification_loss": loss_output.extra_metrics["classification_loss"],
            }
            loss_output_objs.append(loss_dict)

        ### Graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)

        ### Graph KL - mu_all and mu_logvar is the same for both modalities
        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]

        # Log individual graph losses
        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch.values())
        self.log(
            "nll_graph",
            graph_likelihood_loss,
            batch_size=total_batch_size,
            on_epoch=True,
            on_step=False,
        )
        self.log(
            "kl_graph",
            graph_kl_loss_norm,
            batch_size=total_batch_size,
            on_epoch=True,
            on_step=False,
        )

        # Graph loss
        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        # Data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        # Classification loss
        classification_loss = sum(i["classification_loss"] for i in loss_output_objs)
        self.log(
            "class_loss", classification_loss, batch_size=batch_size, on_epoch=True, on_step=False
        )
        if classification_loss > 0 and self.lam_class is None:
            self.lam_class = 100
        elif classification_loss == 0 and self.lam_class is None:
            self.lam_class = 0

        # UOT loss (with annealing if specified)
        z1 = loss_output_objs[0]["z"]
        z2 = loss_output_objs[1]["z"]
        max_epochs = self.trainer.max_epochs
        if self.loss_annealing:
            blur = _anneal_param(
                self.current_epoch, max_epochs, self.init_blur, self.sinkhorn_blur
            )
            reach = _anneal_param(
                self.current_epoch, max_epochs, self.init_reach, self.sinkhorn_reach
            )
        else:
            blur = self.sinkhorn_blur
            reach = self.sinkhorn_reach
        sinkhorn = geomloss.SamplesLoss(loss="sinkhorn", p=self.sinkhorn_p, blur=blur, reach=reach)
        sinkhorn_loss = sinkhorn(z1, z2)
        self.log("uot_loss", sinkhorn_loss, batch_size=batch_size, on_epoch=True, on_step=False)

        # Total loss (lam_kl and lam_data are already included in data_loss)
        total_loss = (
            self.lam_graph * graph_loss
            + data_loss
            + self.lam_sinkhorn * sinkhorn_loss
            + self.lam_class * classification_loss
        )
        self.log(
            "training_loss",
            total_loss,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=total_batch_size,
        )

        return {"loss": total_loss}

    @dependencies("geomloss")
    def validation_step(self, batch: dict[str, dict[str, torch.Tensor]]) -> None:
        """Validation step.
        
        During validation, computes and logs the losses for each modality (NLL, KL, and classification), 
        the graph loss, and the Sinkhorn loss between modalities, and combines them into a total loss.

        Parameters
        ----------
        batch
            A batch of data containing tensors for each modality.

        Returns
        -------
        None
            Returns none. Logs validation losses.
        """
        import geomloss

        # Compute losses for each modality
        loss_output_objs = []
        for _i, (name, tensors) in enumerate(batch.items()):
            batch_size = tensors[REGISTRY_KEYS.X_KEY].shape[0]

            # Update loss kwargs for current modality
            self.loss_kwargs.update(
                {"lam_kl": self.lam_kl, "lam_data": self.lam_data, "mode": name}
            )
            inference_kwargs = {"mode": name}
            generative_kwargs = {"mode": name}

            # Calculate reconstruction, KL, and classification losses
            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            loss = loss_output.loss
            self.log(
                f"val_loss_{name}",
                loss,
                batch_size=batch_size,
                on_epoch=True,
                on_step=True,
            )

            loss_dict = {
                "z": loss_output.extra_metrics["z"],
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
                "classification_loss": loss_output.extra_metrics["classification_loss"],
            }
            loss_output_objs.append(loss_dict)

        # Graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)

        # Graph KL - mu_all and mu_logvar is the same for both modalities
        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]

        # Graph loss
        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        # Data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        # Classification loss
        classification_loss = sum(i["classification_loss"] for i in loss_output_objs)

        # UOT loss
        z1 = loss_output_objs[0]["z"]
        z2 = loss_output_objs[1]["z"]
        sinkhorn = geomloss.SamplesLoss(
            loss="sinkhorn", p=self.sinkhorn_p, blur=self.sinkhorn_blur, reach=self.sinkhorn_reach
        )
        sinkhorn_loss = sinkhorn(z1, z2)

        # Total loss (lam_kl and lam_data are already included in data_loss)
        total_loss = (
            self.lam_graph * graph_loss
            + data_loss
            + self.lam_sinkhorn * sinkhorn_loss
            + self.lam_class * classification_loss
        )
        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch.values())
        self.log(
            "validation_loss",
            total_loss,
            prog_bar=True,
            on_epoch=True,
            batch_size=total_batch_size,
        )
