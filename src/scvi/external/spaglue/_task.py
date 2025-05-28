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

    # edge_weight = graph.edge_weight
    # edge_index = graph.edge_index
    # edge_sign = graph.edge_sign
    # num_nodes = feature_embeddings.size(0)

    # probs = edge_weight / edge_weight.sum()
    # num_samples = edge_index.size(1)
    # # randomly samples from [0, len(probs)-1]
    # sampled_indices = torch.multinomial(probs, num_samples, replacement=True)

    # pos_i = edge_index[0, sampled_indices]
    # pos_j = edge_index[1, sampled_indices]
    # pos_s = edge_sign[sampled_indices]

    # vi = feature_embeddings[pos_i]
    # vj = feature_embeddings[pos_j]

    # pos_logits = pos_s * (vi * vj).sum(dim=1)
    # pos_loss = F.logsigmoid(pos_logits).mean()

    # num_negative = 1

    # # Negative sampling
    # # draw num_samples random ints between 0 and num nodes
    # neg_j = torch.randint(0, num_nodes, (num_negative * num_samples,))
    # neg_i = pos_i.repeat(num_negative)  # if num_negative = 1 just copy the original vector
    # neg_s = pos_s.repeat(num_negative)

    # vi_neg = feature_embeddings[neg_i]
    # vj_neg = feature_embeddings[neg_j]

    # neg_logits = neg_s * (vi_neg * vj_neg).sum(dim=1)
    # neg_loss = F.logsigmoid(-neg_logits).mean()

    # total_loss = -(pos_loss + neg_loss) / 2  # negate the likelihood to get the loss!
    # compute avg
    return total_loss


def kl_divergence_graph(mu, logvar):
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)  # sum over latent dims
    kl_mean = kl.mean()
    return kl_mean


class SPAGLUETrainingPlan(TrainingPlan):
    def __init__(self, *args: tuple, **kwargs: dict) -> None:
        super().__init__(*args, **kwargs)

        self.train_losses = []
        self.validation_losses = []

        self.automatic_optimization = False  # important for adversarial setup

    def training_step(self, batch: list[dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
        # batch contains output of both dataloaders
        # There is a possibility to give batch_idx argument
        # (can be used e.g. for gradient accumulation)
        """Training step."""
        opt = self.optimizers()

        # batch contains both data loader outputs (list of 2 train_dl)
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

            if i == 0:
                self.logger.experiment.add_scalars(
                    "nll", {"seq": reconstruction_loss}, self.global_step
                )
            elif i == 1:
                self.logger.experiment.add_scalars(
                    "nll", {"spatial": reconstruction_loss}, self.global_step
                )

            kl_divergence = loss_output.kl_local["kl_local"]
            kl_divergence = torch.sum(kl_divergence) / (n_obs * n_var)

            if i == 0:
                self.logger.experiment.add_scalars("kl", {"seq": kl_divergence}, self.global_step)
            elif i == 1:
                self.logger.experiment.add_scalars(
                    "kl", {"spatial": kl_divergence}, self.global_step
                )

            loss = loss_output.loss

            loss_dict = {
                "modality_reconstruction_loss": reconstruction_loss,
                "modality_kl_divergence": kl_divergence,
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
                "graph_v_mu": loss_output.extra_metrics["mu_all"],
                "graph_v_logvar": loss_output.extra_metrics["logvar_all"],
            }

            loss_output_objs.append(loss_dict)

        assert torch.equal(loss_output_objs[0]["graph_v_mu"], loss_output_objs[1]["graph_v_mu"])

        ### graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)
        self.logger.experiment.add_scalars(
            "nll", {"graph": graph_likelihood_loss}, self.global_step
        )

        ### graph kl
        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        # in glue: extra normalization with batch size
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]
        self.logger.experiment.add_scalars("kl", {"graph": graph_kl_loss_norm}, self.global_step)

        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        ### data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        K = 1
        total_loss = K * graph_loss + data_loss

        self.train_losses.append(total_loss.item())

        opt.zero_grad()  # zero out gradient
        self.manual_backward(total_loss)  # recompute gradients via backpropagation
        opt.step()  # updates the model weights

        return {"loss": total_loss}

    def validation_step(self, batch: list[dict[str, torch.Tensor]]) -> None:
        """Validation step."""
        print("validating")
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

            # per modality nll
            reconstruction_loss = loss_output.reconstruction_loss["reconstruction_loss"]
            reconstruction_loss = torch.mean(reconstruction_loss)

            # per modality kl
            kl_divergence = loss_output.kl_local["kl_local"]
            kl_divergence = torch.sum(kl_divergence) / (n_obs * n_var)

            loss = loss_output.loss

            loss_dict = {
                "modality_reconstruction_loss": reconstruction_loss,
                "modality_kl_divergence": kl_divergence,
                "modality_loss": loss,
                "graph_v": loss_output.extra_metrics["v_all"],
                "graph_v_mu": loss_output.extra_metrics["mu_all"],
                "graph_v_logvar": loss_output.extra_metrics["logvar_all"],
            }

            loss_output_objs.append(loss_dict)

        ### graph nll
        graph = loss_output.extra_metrics["guidance_graph"]
        feature_embeddings = loss_output_objs[0]["graph_v"]
        graph_likelihood_loss = compute_graph_loss(graph, feature_embeddings)

        ### graph kl
        graph_kl_loss = kl_divergence_graph(
            loss_output.extra_metrics["mu_all"],
            loss_output.extra_metrics["logvar_all"],
        )
        # in glue: extra normalization with batch size
        graph_kl_loss_norm = graph_kl_loss / feature_embeddings.shape[0]

        graph_loss = graph_likelihood_loss + graph_kl_loss_norm

        ### data loss
        data_loss = sum(i["modality_loss"] for i in loss_output_objs)

        K = 1
        total_loss = K * graph_loss + data_loss

        self.validation_losses.append(total_loss.item())
        self.log("validation_loss", total_loss, prog_bar=True, logger=False, on_epoch=True)

    def on_train_epoch_end(self) -> None:
        """Log training loss at the end of each epoch."""
        avg_train_loss = sum(self.train_losses) / len(self.train_losses)
        self.logger.experiment.add_scalars("loss", {"train": avg_train_loss}, self.current_epoch)
        self.train_losses.clear()

    def on_validation_epoch_end(self) -> None:
        """Log validation loss at the end of each epoch."""
        avg_validation_loss = sum(self.validation_losses) / len(self.train_losses)
        self.logger.experiment.add_scalars(
            "loss", {"validation": avg_validation_loss}, self.current_epoch
        )
        self.validation_losses.clear()
