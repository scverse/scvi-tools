"""Training plan for DIAGVI model."""

from __future__ import annotations

from typing import Literal

import torch
import logging

from scvi import REGISTRY_KEYS
from scvi.external.diagvi._utils import (
    compute_graph_loss,
    kl_divergence_graph,
)
from scvi.train import TrainingPlan
from scvi.utils import dependencies

try:
    from pykeops.torch import LazyTensor

    keops_available = True
except ImportError:
    keops_available = False

logger = logging.getLogger(__name__)


def squared_distances(x: torch.Tensor, y: torch.Tensor, use_keops: bool = False) -> torch.Tensor:
    """Compute squared Euclidean distances between point sets."""
    if use_keops and keops_available:
        if x.dim() == 2:
            x_i = LazyTensor(x[:, None, :])
            y_j = LazyTensor(y[None, :, :])
        elif x.dim() == 3:
            x_i = LazyTensor(x[:, :, None, :])
            y_j = LazyTensor(y[:, None, :, :])
        else:
            raise ValueError(f"Incorrect number of dimensions: {x.shape}")
        return ((x_i - y_j) ** 2).sum(-1)

    if x.dim() == 2:
        D_xx = (x * x).sum(-1).unsqueeze(1)
        D_xy = torch.matmul(x, y.permute(1, 0))
        D_yy = (y * y).sum(-1).unsqueeze(0)
    elif x.dim() == 3:
        D_xx = (x * x).sum(-1).unsqueeze(2)
        D_xy = torch.matmul(x, y.permute(0, 2, 1))
        D_yy = (y * y).sum(-1).unsqueeze(1)
    else:
        raise ValueError(f"Incorrect number of dimensions: {x.shape}")

    return D_xx - 2 * D_xy + D_yy


def distances(x: torch.Tensor, y: torch.Tensor, use_keops: bool = False) -> torch.Tensor:
    """Compute Euclidean distances between point sets."""
    if use_keops and keops_available:
        return squared_distances(x, y, use_keops=use_keops).sqrt()
    return torch.sqrt(torch.clamp_min(squared_distances(x, y), 1e-8))


cost_routines = {
    1: lambda x, y: distances(x, y),
    2: lambda x, y: squared_distances(x, y) / 2,
}


def _compute_cost_scale(
    C: torch.Tensor,
    method: float | Literal["mean", "max", "median", "std"] | None,
) -> float:
    """Compute scale factor for cost matrix (OTT-JAX style).

    Returns the scale factor s such that C_scaled = C / s has unit statistic.
    """
    if method is None or method == 1.0:
        return 1.0
    if isinstance(method, (int | float)):
        return float(method)
    if method == "mean":
        return C.mean().item()
    if method == "max":
        return C.max().item()
    if method == "median":
        return C.median().item()
    if method == "std":
        return C.std().item()
    raise ValueError(f"Unknown scale_cost method: {method}")


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
    Annealed parameter value.
    """
    anneal_epochs = max_epochs // 3
    if current_epoch >= anneal_epochs:
        return target_value
    progress = current_epoch / anneal_epochs
    return init_value + (target_value - init_value) * progress


class DiagTrainingPlan(TrainingPlan):
    """Training plan for DIAGVI model.

    Handles multi-modal loss computation including reconstruction, KL divergence,
    graph loss, Sinkhorn (optimal transport) loss, and classification loss.

    Parameters
    ----------
    module
        DIAGVAE module to train.
    lam_graph
        Weight for the graph reconstruction loss.
    lam_kl
        Weight for the KL divergence loss.
    lam_data
        Weight for the data reconstruction loss.
    lam_sinkhorn
        Weight for the Sinkhorn (optimal transport) loss.
    lam_class
        Weight for the classification loss. If None, automatically set based
        on whether classification is used.
    sinkhorn_p
        Order of the Wasserstein distance (p-norm).
    sinkhorn_blur
        Blur for geomloss Sinkhorn (epsilon = blur^p). If None, computed
        adaptively from cost matrix like OTT-JAX, by default None.
    epsilon_from_cost
        Statistic of scaled cost to compute epsilon: "std" or "mean", by default "mean".
    epsilon_scale
        Scaling factor: `epsilon = epsilon_scale * statistic(C_scaled)`, by default 0.05.
    scale_cost
        How to scale cost matrix. Either a float (C_scaled = C / float) or
        one of "mean", "max", "median", "std" (C_scaled = C / statistic(C)).
        Use None or 1.0 for no scaling. By default 1.0.
    sinkhorn_reach
        Reach parameter for unbalanced OT. If None, calculate adaptive reach from blur.
    reach_scale
        Scaling factor for adaptive reach: `reach = reach_scale * blur`, by default 10.0.
    lr
        Learning rate.
    loss_annealing
        Whether to anneal loss parameters during training.
    log_train
        Whether to log individual train loss components.
    log_val
        Whether to log individual validation loss components.
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
        lam_class: float = 100.0,
        sinkhorn_p: int = 2,
        sinkhorn_blur: float | None = None,
        epsilon_from_cost: Literal["mean", "std"] = "mean",
        epsilon_scale: float = 0.05,
        scale_cost: float | Literal["mean", "max", "median", "std"] | None = 1.0,
        sinkhorn_reach: float | None = None,
        reach_scale: float = 10.0,
        lr: float = 1e-3,
        loss_annealing: bool = False,
        log_train: bool = True,
        log_val: bool = False,
        *args,
        **kwargs,
    ) -> None:
        super().__init__(module, *args, **kwargs)

        # Loss weights
        self.lam_graph = lam_graph
        self.lam_kl = lam_kl
        self.lam_data = lam_data
        self.lam_sinkhorn = lam_sinkhorn
        self.lam_class = lam_class

        # Sinkhorn parameters
        self.sinkhorn_p = sinkhorn_p
        self.sinkhorn_blur = sinkhorn_blur
        self.sinkhorn_reach = sinkhorn_reach
        self.epsilon_from_cost = epsilon_from_cost
        self.epsilon_scale = epsilon_scale
        self.scale_cost = scale_cost
        self.reach_scale = reach_scale

        # Training parameters
        self.lr = lr
        self.loss_annealing = loss_annealing
        self.log_train = log_train
        self.log_val = log_val

        # Initial values for annealing (10x larger for smoother optimization start)
        self.init_blur = 10 * self.sinkhorn_blur if self.sinkhorn_blur is not None else None
        self.init_reach = 10 * self.sinkhorn_reach if self.sinkhorn_reach is not None else None

    def _compute_modality_losses(
        self, batch: dict[str, dict[str, torch.Tensor]], log_prefix: str = ""
    ) -> tuple[list[dict], dict, int]:
        """Compute per-modality losses (reconstruction + KL).

        Parameters
        ----------
        batch
            A batch of data containing tensors for each modality.
        log_prefix
            Prefix for logging keys (e.g., "train_" or "val_").

        Returns
        -------
        loss_outputs
            Per-modality loss information containing z, modality_loss, graph_v, classification_loss
        last_loss_output
            The last loss output object (contains shared info like guidance_graph)
        total_batch_size
            Total number of samples across all modalities
        """
        loss_outputs = []

        for name, tensors in batch.items():
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
                inference_kwargs={"mode": name},
                generative_kwargs={"mode": name},
            )

            # Log per-modality metrics
            # Validation and training step
            if log_prefix:
                self.log(
                    f"{log_prefix}loss_{name}",
                    loss_output.loss,
                    batch_size=batch_size,
                    on_epoch=True,
                    on_step=True,
                )

            # Only training step
            if log_prefix == "train_":
                reconstruction_loss = torch.mean(
                    loss_output.reconstruction_loss["reconstruction_loss"]
                )
                self.log(f"nll_{name}", reconstruction_loss, batch_size=batch_size, on_epoch=True)
                self.log(
                    f"kl_{name}",
                    loss_output.kl_local["kl_local"],
                    batch_size=batch_size,
                    on_epoch=True,
                )

            loss_outputs.append(
                {
                    "z": loss_output.extra_metrics["z"],
                    "modality_loss": loss_output.loss,
                    "graph_v": loss_output.extra_metrics["v_all"],
                    "classification_loss": loss_output.extra_metrics["classification_loss"],
                }
            )

        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch.values())

        return loss_outputs, loss_output, total_batch_size

    def _compute_graph_loss(
        self,
        loss_output: dict,
        feature_embeddings: torch.Tensor,
        log: bool = False,
        batch_size: int = 0
    ) -> torch.Tensor:
        """Compute graph reconstruction and KL loss.
        
        Parameters
        ----------
        loss_output
            Loss output object containing graph information.
        feature_embeddings
            Tensor of shape (n_features, embedding_dim) containing feature embeddings.
        log
            Whether to log graph loss components.
        batch_size
            Batch size for logging.

        Returns
        -------
        Scalar tensor containing the total graph loss.
        """
        graph = loss_output.extra_metrics["guidance_graph"]

        # Graph reconstruction loss
        graph_nll = compute_graph_loss(graph, feature_embeddings)

        # Graph KL divergence (normalized by number of features)
        graph_kl = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        graph_kl_norm = graph_kl / feature_embeddings.shape[0]

        if log:
            self.log("nll_graph", graph_nll, batch_size=batch_size, on_epoch=True)
            self.log("kl_graph", graph_kl_norm, batch_size=batch_size, on_epoch=True)

        return graph_nll + graph_kl_norm

    @dependencies("geomloss")
    def _compute_sinkhorn_loss(
        self,
        z1: torch.Tensor,
        z2: torch.Tensor,
        use_annealing: bool = False
    ) -> torch.Tensor:
        """Compute Sinkhorn (Unbalanced Optimal Transport) loss between latent spaces.
        
        Parameters
        ----------
        z1
            Latent representations from modality 1, shape (n_cells, n_latent).
        z2
            Latent representations from modality 2, shape (n_cells, n_latent).
        use_annealing
            Whether to use annealed Sinkhorn parameters.

        Returns
        -------
        Scalar tensor containing the Sinkhorn loss.
        """
        import geomloss

        logger.info(
            f"Computing UOT with p={self.sinkhorn_p}, blur={self.sinkhorn_blur}, reach={self.sinkhorn_reach}"
        )
        logger.info(f"Cost scaling method: {self.scale_cost}")
        if self.sinkhorn_blur is None:
            logger.info(
                f"blur not provided, will compute adaptively from cost matrix using epsilon_from_cost={self.epsilon_from_cost} and epsilon_scale={self.epsilon_scale}"
            )
        if self.sinkhorn_reach is None:
            logger.info(
                f"reach not provided, will compute adaptively from blur using reach_scale={self.reach_scale}."
            )

        # Step 1: Compute cost matrix C
        cost_fn = cost_routines[self.sinkhorn_p]
        C_st = cost_fn(z1, z2.detach())

        # Step 2: Compute cost scale factor (OTT-JAX style)
        cost_scale = _compute_cost_scale(C_st, self.scale_cost)
        if cost_scale <= 0:
            cost_scale = 1.0  # Fallback for edge cases

        # Step 3: Scale cost matrix: C_scaled = C / cost_scale
        C_scaled = C_st / cost_scale

        # Step 4: Compute adaptive epsilon from SCALED cost (OTT-JAX approach)
        if self.sinkhorn_blur is None:
            if self.epsilon_from_cost == "std":
                eps = self.epsilon_scale * C_scaled.std().item()
            elif self.epsilon_from_cost == "mean":
                eps = self.epsilon_scale * C_scaled.mean().item()
            else:
                raise ValueError(f"Unknown epsilon_from_cost: {self.epsilon_from_cost}")

            eps = max(eps, 1e-8)  # Ensure positive

            # Convert epsilon to blur for geomloss: epsilon = blur^p → blur = epsilon^(1/p)
            blur = eps ** (1.0 / self.sinkhorn_p)
        else:
            eps = self.sinkhorn_blur**self.sinkhorn_p
            blur = self.sinkhorn_blur

        # Step 5: Compute adaptive reach if needed
        if self.sinkhorn_reach is None:
            reach = self.reach_scale * blur  # Heuristic: reach should be larger than blur

        # Step 6: Scale embeddings so geomloss sees scaled cost
        # We want: cost(x_scaled, y_scaled) = C_st / cost_scale = C_scaled
        # For cost = ||x-y||^p / p:
        #   ||x_scaled - y_scaled||^p / p = ||x - y||^p / (p * cost_scale)
        # If x_scaled = x * k, then k^p = 1 / cost_scale → k = cost_scale^(-1/p)
        if cost_scale != 1.0:
            embed_scale = cost_scale ** (-1.0 / self.sinkhorn_p)
            z1_scaled = z1 * embed_scale
            z2_scaled = z2 * embed_scale
        else:
            z1_scaled = z1
            z2_scaled = z2
        
        if use_annealing and self.loss_annealing:
            logger.info("Using annealed Sinkhorn parameters.")
            max_epochs = self.trainer.max_epochs
            blur = _anneal_param(
                self.current_epoch, max_epochs, self.init_blur, self.sinkhorn_blur
            )
            reach = _anneal_param(
                self.current_epoch, max_epochs, self.init_reach, self.sinkhorn_reach
            )
        else:
            blur = blur
            reach = reach
        
        sinkhorn = geomloss.SamplesLoss(loss="sinkhorn", p=self.sinkhorn_p, blur=blur, reach=reach)
        return sinkhorn(z1_scaled, z2_scaled)

    def _compute_total_loss(
        self,
        loss_outputs: list[dict],
        loss_output: dict,
        total_batch_size: int,
        use_annealing: bool = False,
        log: bool = False,
    ) -> torch.Tensor:
        """Compute total loss combining all components.
        
        Parameters
        ----------
        loss_outputs
            List of per-modality loss information.
        loss_output
            The last loss output object (contains shared info like guidance_graph).
        total_batch_size
            Total number of samples across all modalities.
        use_annealing
            Whether to use annealed Sinkhorn parameters.
        log
            Whether to log loss components.

        Returns
        -------
        Scalar tensor containing the total loss.
        """
        # 1. Data loss (reconstruction + KL for each modality)
        data_loss = sum(out["modality_loss"] for out in loss_outputs)

        # 2. Graph loss
        feature_embeddings = loss_outputs[0]["graph_v"]
        graph_loss = self._compute_graph_loss(
            loss_output, feature_embeddings, log=log, batch_size=total_batch_size
        )

        # 3. Classification loss
        classification_loss = sum(out["classification_loss"] for out in loss_outputs)
        if log:
            self.log("class_loss", classification_loss, batch_size=total_batch_size, on_epoch=True)

        # 4. Sinkhorn (UOT) loss
        z1, z2 = loss_outputs[0]["z"], loss_outputs[1]["z"]
        sinkhorn_loss = self._compute_sinkhorn_loss(z1, z2, use_annealing=use_annealing)
        if log:
            self.log("uot_loss", sinkhorn_loss, batch_size=total_batch_size, on_epoch=True)

        # Combine all losses
        total_loss = (
            data_loss
            + self.lam_graph * graph_loss
            + self.lam_sinkhorn * sinkhorn_loss
            + self.lam_class * classification_loss
        )

        return total_loss

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
        A dictionary containing the total loss for backpropagation.
        """
        loss_outputs, loss_output, total_batch_size = self._compute_modality_losses(
            batch, log_prefix="train_"
        )

        total_loss = self._compute_total_loss(
            loss_outputs, loss_output, total_batch_size, use_annealing=self.loss_annealing, log=self.log_train
        )

        self.log(
            "training_loss",
            total_loss,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=total_batch_size,
        )

        # Return embeddings along with loss for callbacks (detached to avoid graph issues)
        embeddings = {
            name: out["z"].detach()
            for name, out in zip(batch.keys(), loss_outputs, strict=False)
        }

        return {"loss": total_loss, "embeddings": embeddings}

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
        Returns none. Logs validation losses.
        """
        loss_outputs, loss_output, total_batch_size = self._compute_modality_losses(
            batch, log_prefix="val_"
        )

        total_loss = self._compute_total_loss(
            loss_outputs, loss_output, total_batch_size, use_annealing=False, log=self.log_val
        )

        self.log(
            "validation_loss",
            total_loss,
            prog_bar=True,
            on_epoch=True,
            batch_size=total_batch_size,
        )
