import torch
import torch.nn.functional as F
from torch_geometric.utils import structured_negative_sampling

from scvi import REGISTRY_KEYS
from scvi.train import TrainingPlan


def compute_graph_loss(graph, feature_embeddings):
    edge_index = graph.edge_index
    edge_index_neg = structured_negative_sampling(edge_index)

    pos_i = edge_index_neg[0].cpu().numpy()
    pos_j = edge_index_neg[1].cpu().numpy()
    neg_j = edge_index_neg[2].cpu().numpy()

    vi = feature_embeddings[pos_i]
    vj = feature_embeddings[pos_j]
    vj_neg = feature_embeddings[neg_j]

    pos_logits = (vi * vj).sum(dim=1)
    pos_loss = F.logsigmoid(pos_logits).mean()

    neg_logits = (vi * vj_neg).sum(dim=1)
    neg_loss = F.logsigmoid(-neg_logits).mean()

    total_loss = -(pos_loss + neg_loss) / 2

    return total_loss


def kl_divergence_graph(mu, logvar):
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)  # sum over latent dims
    kl_mean = kl.mean()
    return kl_mean


class SPAGLUETrainingPlan(TrainingPlan):
    def __init__(self, module, lam_graph=1.0, lam_kl=1.0, *args, **kwargs) -> None:
        super().__init__(module, *args, **kwargs)

        self.lam_graph = lam_graph
        self.lam_kl = lam_kl

        self.automatic_optimization = False  # important for adversarial setup

    def training_step(self, batch: list[dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
        # batch contains output of both dataloaders
        # There is a possibility to give batch_idx argument
        # (can be used e.g. for gradient accumulation)
        """Training step."""
        opt = self.optimizers()
        loss_output_objs = []

        for i, tensors in enumerate(batch):
            batch_size = tensors[REGISTRY_KEYS.X_KEY].shape[0]

            self.loss_kwargs.update({"lam_kl": self.lam_kl, "mode": i})
            inference_kwargs = {"mode": i}
            generative_kwargs = {"mode": i}

            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            # just for logging
            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)
            if i == 0:
                self.log("nll_seq", reconstruction_loss, batch_size=batch_size)
            elif i == 1:
                self.log("nll_spatial", reconstruction_loss, batch_size=batch_size)
            kl_divergence = loss_output.kl_local["kl_local"]
            if i == 0:
                self.log("kl_seq", kl_divergence, batch_size=batch_size)
            elif i == 1:
                self.log("kl_spatial", kl_divergence, batch_size=batch_size)

            loss = loss_output.loss

            loss_dict = {
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
            }

            loss_output_objs.append(loss_dict)

        ### graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)

        ### graph kl - mu_all and mu_logvar is the same for both modalities
        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]

        # log individual graph losses
        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch)
        self.log("nll_graph", graph_likelihood_loss, batch_size=total_batch_size)
        self.log("kl_graph", graph_kl_loss_norm, batch_size=total_batch_size)

        ### graph loss
        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        ### data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        ### total loss
        total_loss = self.lam_graph * graph_loss + data_loss

        self.log(
            "training_loss",
            total_loss,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=total_batch_size,
        )

        opt.zero_grad()
        self.manual_backward(total_loss)
        opt.step()
        return {"loss": total_loss}

    def validation_step(self, batch: list[dict[str, torch.Tensor]]) -> None:
        """Validation step."""
        loss_output_objs = []

        for i, tensors in enumerate(batch):
            self.loss_kwargs.update({"lam_kl": self.lam_kl, "mode": i})
            inference_kwargs = {"mode": i}
            generative_kwargs = {"mode": i}

            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            loss = loss_output.loss

            loss_dict = {
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
            }

            loss_output_objs.append(loss_dict)

            ### graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]  # 0 or 1 is same
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)

        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]

        ### graph loss
        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        ### data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        ### total loss
        total_loss = self.lam_graph * graph_loss + data_loss

        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch)
        self.log(
            "validation_loss",
            total_loss,
            prog_bar=True,
            on_epoch=True,
            batch_size=total_batch_size,
        )
