import geomloss
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


def distance_matrix(pts_src: torch.Tensor, pts_dst: torch.Tensor, p: int = 2):
    x_col = pts_src.unsqueeze(1)
    y_row = pts_dst.unsqueeze(0)
    distance = torch.sum((torch.abs(x_col - y_row)) ** p, 2)
    return distance


def kl_divergence_graph(mu, logvar):
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)  # sum over latent dims
    kl_mean = kl.mean()
    return kl_mean


class SPAGLUETrainingPlan(TrainingPlan):
    def __init__(
        self,
        module,
        lam_graph=1.0,
        lam_kl=1.0,
        lam_data=1.0,
        lam_sinkhorn=1.0,
        sinkhorn_p=2,
        sinkhorn_blur=1,
        sinkhorn_reach=1,
        lr=1e-3,
        *args,
        **kwargs,
    ) -> None:
        super().__init__(module, *args, **kwargs)

        self.lam_graph = lam_graph
        self.lam_kl = lam_kl
        self.lam_data = lam_data
        self.lam_sinkhorn = lam_sinkhorn
        self.sinkhorn_p = sinkhorn_p
        self.sinkhorn_reach = sinkhorn_reach
        self.sinkhorn_blur = sinkhorn_blur
        self.lr = lr  # scvi handles giving the learning rate to the optimizer

        # self.automatic_optimization = False

    def training_step(self, batch: dict[str, dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
        """Training step."""
        # opt = self.optimizers()
        # print(opt)

        loss_output_objs = []
        for _i, (modality, tensors) in enumerate(batch.items()):
            batch_size = tensors[REGISTRY_KEYS.X_KEY].shape[0]

            self.loss_kwargs.update(
                {"lam_kl": self.lam_kl, "lam_data": self.lam_data, "mode": modality}
            )
            inference_kwargs = {"mode": modality}
            generative_kwargs = {"mode": modality}

            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            # just for logging
            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)
            self.log("nll_{modality}", reconstruction_loss, batch_size=batch_size)

            kl_divergence = loss_output.kl_local["kl_local"]
            self.log("kl_{modality}", kl_divergence, batch_size=batch_size)

            loss = loss_output.loss

            loss_dict = {
                "z": loss_output.extra_metrics["z"],
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
        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch.values())
        self.log("nll_graph", graph_likelihood_loss, batch_size=total_batch_size)
        self.log("kl_graph", graph_kl_loss_norm, batch_size=total_batch_size)

        ### graph loss
        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        ### data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        ### UOT loss
        z1 = loss_output_objs[0]["z"]
        z2 = loss_output_objs[1]["z"]

        # uot_loss, tran = unbalanced_ot(z1, z2)
        sinkhorn = geomloss.SamplesLoss(
            loss="sinkhorn", p=self.sinkhorn_p, blur=self.sinkhorn_blur, reach=self.sinkhorn_reach
        )
        sinkhorn_loss = sinkhorn(z1, z2)

        self.log("uot_loss", sinkhorn_loss, batch_size=batch_size)

        ### total loss
        total_loss = self.lam_graph * graph_loss + data_loss + self.lam_sinkhorn * sinkhorn_loss

        self.log(
            "training_loss",
            total_loss,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=total_batch_size,
        )

        # opt.zero_grad()
        # self.manual_backward(total_loss)
        # opt.step()

        return {"loss": total_loss}

    def validation_step(self, batch: list[dict[str, torch.Tensor]]) -> None:
        """Validation step."""
        loss_output_objs = []

        for _i, (modality, tensors) in enumerate(batch.items()):
            self.loss_kwargs.update(
                {"lam_kl": self.lam_kl, "lam_data": self.lam_data, "mode": modality}
            )
            inference_kwargs = {"mode": modality}
            generative_kwargs = {"mode": modality}

            _, _, loss_output = self.forward(
                tensors,
                loss_kwargs=self.loss_kwargs,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
            )

            loss = loss_output.loss

            loss_dict = {
                "z": loss_output.extra_metrics["z"],
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

        ### UOT loss
        z1 = loss_output_objs[0]["z"]
        z2 = loss_output_objs[1]["z"]

        # uot_loss, tran = unbalanced_ot(z1, z2)
        sinkhorn = geomloss.SamplesLoss(
            loss="sinkhorn", p=self.sinkhorn_p, blur=self.sinkhorn_blur, reach=self.sinkhorn_reach
        )
        sinkhorn_loss = sinkhorn(z1, z2)

        ### total loss (lam_kl and lam_data are already included in data_loss)
        total_loss = self.lam_graph * graph_loss + data_loss + self.lam_sinkhorn * sinkhorn_loss

        total_batch_size = sum(tensors[REGISTRY_KEYS.X_KEY].shape[0] for tensors in batch.values())
        self.log(
            "validation_loss",
            total_loss,
            prog_bar=True,
            on_epoch=True,
            batch_size=total_batch_size,
        )
