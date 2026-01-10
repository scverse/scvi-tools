import logging

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical, Independent, MixtureSameFamily, kl_divergence

from scvi import REGISTRY_KEYS

logger = logging.getLogger(__name__)
from scvi.distributions import (
    Gamma,
    Normal,
    LogNormal,
    Log1pNormal,
    NegativeBinomial,
    NegativeBinomialMixture,
    ZeroInflatedGamma,
    ZeroInflatedLogNormal,
    ZeroInflatedNegativeBinomial,
)
from scvi.external.diagvi import DecoderDualPathway, DecoderSinglePathway, GraphEncoder_glue
from scvi.module import Classifier
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder

# Decoder selection based on likelihood function
# Single-pathway: outputs (scale, r, rate, dropout)
# Dual-pathway: outputs ((scale1, scale2), r, (rate1, rate2), mixture_logits)
LIKELIHOOD_TO_DECODER = {
    # Single-pathway likelihoods (unimodal distributions)
    "nb": DecoderSinglePathway,
    "zinb": DecoderSinglePathway,
    "normal": DecoderSinglePathway,
    "lognormal": DecoderSinglePathway,
    "log1pnormal": DecoderSinglePathway,
    "ziln": DecoderSinglePathway,
    "gamma": DecoderSinglePathway,
    "zig": DecoderSinglePathway,
    # Dual-pathway likelihoods (mixture distributions)
    "nbmixture": DecoderDualPathway,
}

# Likelihoods that require softmax normalization (count data)
COUNT_LIKELIHOODS = {"nb", "zinb", "nbmixture"}


class DIAGVAE(BaseModuleClass):
    """
    Variational autoencoder module for DIAGVI multi-modal integration.

    Supports GMM priors, semi-supervised classification, and flexible modality-specific decoders.

    Parameters
    ----------
    n_inputs : dict[str, int]
        Number of input features for each modality.
    n_batches : dict[str, int]
        Number of batches for each modality.
    n_labels : dict[str, int]
        Number of labels/classes for each modality.
    modality_likelihoods : dict[str, str]
        Likelihood model for each modality (e.g., 'nb', 'zinb', 'nbmixture', 'normal').
    modalities : dict[str, str]
        Modality type for each input (e.g., 'rna', 'protein').
    guidance_graph : torch_geometric.data.Data
        Graph object encoding feature correspondences.
    use_gmm_prior : dict[str, bool]
        Whether to use a GMM prior for each modality.
    semi_supervised : dict[str, bool]
        Whether to use semi-supervised classification for each modality.
    n_mixture_components : dict[str, int]
        Number of mixture components for the GMM prior for each modality.
    n_latent : int, optional
        Dimensionality of the latent space (default: 50).
    n_hidden : int, optional
        Number of nodes per hidden layer (default: 256).
    n_layers : int, optional
        Number of hidden layers (default: 2).
    dropout_rate : float, optional
        Dropout rate for encoders (default: 0.1).
    """

    def __init__(
        self,
        n_inputs: dict[int],
        n_batches: dict[int],
        n_labels: dict[int],
        modality_likelihoods: dict[str],
        modalities: dict[str],
        guidance_graph,
        use_gmm_prior: dict[bool],
        semi_supervised: dict[bool],
        n_mixture_components: dict[int],
        n_latent: int = 50,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        common_scale: bool = True,
        # **kwargs: dict,
    ) -> None:
        super().__init__()

        # learnable parameters
        self.gmm_logits = nn.ParameterDict()
        self.gmm_means = nn.ParameterDict()
        self.gmm_scales = nn.ParameterDict()

        self.use_gmm_prior = use_gmm_prior
        self.semi_supervised = semi_supervised
        self.n_mixture_components = n_mixture_components
        self.n_labels = n_labels
        self.n_input_list = n_inputs
        self.n_batches_list = n_batches
        self.modality_likelihoods = modality_likelihoods
        self.modalities = modalities
        self.guidance_graph = guidance_graph
        self.n_latent = n_latent

        self.input_names = list(n_inputs.keys())

        for name in use_gmm_prior.keys():
            if self.use_gmm_prior[name]:
                k = self.n_mixture_components[name]
                # overwrite the number of mixture components if semi_supervised = True
                if self.semi_supervised[name]:
                    k = self.n_labels[name]
                self.gmm_logits[name] = nn.Parameter(torch.zeros(k))
                self.gmm_means[name] = nn.Parameter(torch.randn(k, n_latent))
                self.gmm_scales[name] = nn.Parameter(torch.zeros(k, n_latent))


        # encoders
        self.encoder_0 = Encoder(
            n_input=n_inputs[self.input_names[0]],
            n_output=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )

        self.encoder_1 = Encoder(
            n_input=n_inputs[self.input_names[1]],
            n_output=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )


        # decoders - selected based on likelihood function
        likelihood_0 = modality_likelihoods[self.input_names[0]]
        likelihood_1 = modality_likelihoods[self.input_names[1]]

        # Determine decoder class based on likelihood
        decoder_class_0 = LIKELIHOOD_TO_DECODER.get(likelihood_0, DecoderSinglePathway)
        decoder_class_1 = LIKELIHOOD_TO_DECODER.get(likelihood_1, DecoderSinglePathway)

        # Determine whether to apply softmax normalization based on likelihood
        normalize_0 = likelihood_0 in COUNT_LIKELIHOODS
        normalize_1 = likelihood_1 in COUNT_LIKELIHOODS

        logger.info(
            f"Decoder for '{self.input_names[0]}' (likelihood={likelihood_0}): "
            f"{decoder_class_0.__name__}, normalize={normalize_0}"
        )
        self.decoder_0 = decoder_class_0(
            n_output=n_inputs[self.input_names[0]],
            n_batches=n_batches[self.input_names[0]],
            normalize=normalize_0,
        )

        logger.info(
            f"Decoder for '{self.input_names[1]}' (likelihood={likelihood_1}): "
            f"{decoder_class_1.__name__}, normalize={normalize_1}"
        )
        self.decoder_1 = decoder_class_1(
            n_output=n_inputs[self.input_names[1]],
            n_batches=n_batches[self.input_names[1]],
            normalize=normalize_1,
        )

        self.graph_encoder = GraphEncoder_glue(
            vnum=n_inputs[self.input_names[0]] + n_inputs[self.input_names[1]],
            out_features=self.n_latent,
        )

        if self.semi_supervised[self.input_names[0]]:
            cls_parameters = {
                "n_layers": 0,
                "n_hidden": 128,
                "dropout_rate": 0.0,
            }
            self.classifier_0 = Classifier(
                n_latent,
                n_labels=self.n_labels[self.input_names[0]],
                use_batch_norm=False,
                use_layer_norm=True,
                **cls_parameters,
            )
        else:
            self.classifier_0 = None

        if self.semi_supervised[self.input_names[1]]:
            cls_parameters = {
                "n_layers": 0,
                "n_hidden": 128,
                "dropout_rate": 0.0,
            }
            self.classifier_1 = Classifier(
                n_latent,
                n_labels=self.n_labels[self.input_names[1]],
                use_batch_norm=False,
                use_layer_norm=True,
                **cls_parameters,
            )
        else:
            self.classifier_1 = None

    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor | None]:
        """Prepare input dictionary for the inference step."""
        return {
            MODULE_KEYS.X_KEY: tensors[REGISTRY_KEYS.X_KEY],
        }

    def _get_generative_input(
        self, tensors: dict[str, torch.Tensor], inference_outputs: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor]:
        """Prepare input dictionary for the generative step."""
        return {
            MODULE_KEYS.Z_KEY: inference_outputs[MODULE_KEYS.Z_KEY],
            MODULE_KEYS.LIBRARY_KEY: inference_outputs[MODULE_KEYS.LIBRARY_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            MODULE_KEYS.Y_KEY: tensors[REGISTRY_KEYS.LABELS_KEY],
            "v": inference_outputs["v"],
        }

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        mode: str | None = None,
    ) -> dict[str, torch.Tensor]:
        """
        Run the inference (encoder and graph) step for a given modality.

        Parameters
        ----------
        x : torch.Tensor
            Input data tensor.
        mode : str, optional
            Name of the modality.

        Returns
        -------
        dict[str, torch.Tensor]
            Dictionary of inference outputs, including latent variables and graph embeddings.
        """
        x_ = x
        library = torch.log(x.sum(1)).unsqueeze(1)
        graph = self.guidance_graph
        device = x.device
        graph = graph.to(device)
        # graph inference
        v_all, mu_all, logvar_all = self.graph_encoder(graph.edge_index)
        v = v_all[getattr(graph, f"{mode}_indices")]
        other_mode = [m for m in self.input_names if m != mode][0]
        v_other_mod = v_all[getattr(graph, f"{other_mode}_indices")]
        # data inference
        if mode == self.input_names[0]:
            qz, z = self.encoder_0(x_)
        else:
            qz, z = self.encoder_1(x_)
        return {
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.LIBRARY_KEY: library,
            "v": v,
            "v_other": v_other_mod,
            "v_all": v_all,
            "mu_all": mu_all,
            "logvar_all": logvar_all,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        y: torch.Tensor | None = None,
        v: torch.Tensor | None = None,
        mode: str | None = 0,
    ) -> dict[str, torch.Tensor]:
        """
        Run the generative (decoder) model for a given modality.

        Parameters
        ----------
        z : torch.Tensor
            Latent variable tensor.
        library : torch.Tensor
            Library size tensor.
        batch_index : torch.Tensor, optional
            Batch index tensor.
        y : torch.Tensor, optional
            Cell type labels (for semi-supervised GMM prior).
        v : torch.Tensor, optional
            Feature embedding from the graph encoder.
        mode : str, optional
            Name of the modality.

        Returns
        -------
        dict[str, torch.Tensor]
            Dictionary containing generative model outputs and distributions.
        """
        EPS = 1e-8
        if mode == self.input_names[0]:
            px_scale, px_r, px_rate, px_dropout = self.decoder_0(z, library, batch_index, v)
        elif mode == self.input_names[1]:
            px_scale, px_r, px_rate, px_dropout = self.decoder_1(z, library, batch_index, v)
        if self.modality_likelihoods[mode] in COUNT_LIKELIHOODS:
            px_r = px_r.exp()
        else:
            px_r = F.softplus(px_r) + EPS

        if self.modality_likelihoods[mode] == "nb":
            px = NegativeBinomial(px_r, logits=(px_rate + EPS).log() - px_r)
        elif self.modality_likelihoods[mode] == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px_rate,
                theta=px_r,
                zi_logits=px_dropout,
                scale=px_scale,
            )
        elif self.modality_likelihoods[mode] == "normal":
            px = Normal(px_rate, px_r, normal_mu=px_scale)
        elif self.modality_likelihoods[mode] == "nbmixture":
            px = NegativeBinomialMixture(
                mu1=px_rate[0],
                mu2=px_rate[1],
                theta1=px_r,
                mixture_logits=px_dropout,
            )
            """
            px = NegativeBinomialMixture(
                mu1=px_scale,
                mu2=px_r,
                theta1=px_rate,
                mixture_logits=px_dropout,
            )
            """
        elif self.modality_likelihoods[mode] == "lognormal":
            px = LogNormal(mu=px_rate, sigma=px_r, scale=px_scale)
        elif self.modality_likelihoods[mode] == "log1pnormal":
            px = Log1pNormal(mu=px_rate, sigma=px_r, scale=px_scale)
        elif self.modality_likelihoods[mode] == "gamma":
            # Clamp concentration and rate to avoid lgamma numerical issues
            px = Gamma(
                concentration=torch.clamp(px_rate, min=EPS),
                rate=torch.clamp(px_r, min=EPS),
            )
        elif self.modality_likelihoods[mode] == "ziln":
            px = ZeroInflatedLogNormal(
                mu=px_rate,
                sigma=px_r,
                zi_logits=px_dropout,
                scale=px_scale,
            )
        elif self.modality_likelihoods[mode] == "zig":
            # Clamp concentration and rate to avoid lgamma numerical issues
            px = ZeroInflatedGamma(
                concentration=torch.clamp(px_rate, min=EPS),
                rate=torch.clamp(px_r, min=EPS),
                zi_logits=px_dropout,
                scale=px_scale,
            )
        else:
            raise ValueError(
                f"Unknown likelihood '{self.modality_likelihoods[mode]}' for modality '{mode}'. "
                f"Supported likelihoods: {list(LIKELIHOOD_TO_DECODER.keys())}"
            )

        if self.use_gmm_prior[mode]:
            logits = self.gmm_logits[mode]
            means = self.gmm_means[mode]
            scales = torch.exp(self.gmm_scales[mode]) + 1e-4
            if self.semi_supervised[mode]:
                logits_input = (
                    torch.stack(
                        [
                            torch.nn.functional.one_hot(y_i, self.n_labels[mode])
                            if y_i < self.n_labels[mode]
                            else torch.zeros(self.n_labels[mode])
                            for y_i in y.ravel()
                        ]
                    )
                    .to(z.device)
                    .float()
                )
                logits = logits + 100 * logits_input
                means = means.expand(y.shape[0], -1, -1)
                scales = scales.expand(y.shape[0], -1, -1)
            cats = Categorical(logits=logits)
            normal_dists = Independent(Normal(means, scales), reinterpreted_batch_ndims=1)
            pz = MixtureSameFamily(cats, normal_dists)
        else:
            pz = Normal(torch.zeros_like(z), torch.ones_like(z))
        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PZ_KEY: pz,
            "px_rate": px_rate,
        }

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        lam_kl: torch.Tensor | float = 1.0,
        lam_data: torch.Tensor | float = 1.0,
        mode: str | None = None,
    ) -> LossOutput:
        """Compute the loss for a batch"""
        x = tensors[REGISTRY_KEYS.X_KEY]
        n_obs = x.shape[0]
        n_var = x.shape[1]

        # data nll calculation
        reconst_loss = -generative_outputs[MODULE_KEYS.PX_KEY].log_prob(x).sum(-1)
        reconstruction_loss_norm = torch.mean(reconst_loss)
        if self.use_gmm_prior[mode]:
            kl_div = inference_outputs[MODULE_KEYS.QZ_KEY].log_prob(inference_outputs["z"]).sum(
                -1
            ) - generative_outputs[MODULE_KEYS.PZ_KEY].log_prob(inference_outputs["z"])
        else:
            kl_div = kl_divergence(
                inference_outputs[MODULE_KEYS.QZ_KEY], generative_outputs[MODULE_KEYS.PZ_KEY]
            ).sum(dim=-1)
        kl_local_norm = torch.sum(kl_div) / (n_obs * n_var)
        loss = lam_data * reconstruction_loss_norm + lam_kl * kl_local_norm

        ## graph inference
        mu_all = inference_outputs["mu_all"]
        logvar_all = inference_outputs["logvar_all"]
        v_all = inference_outputs["v_all"]

        classification_loss = 0.0
        if self.classifier_0 is not None and mode == self.input_names[0]:
            y = tensors[REGISTRY_KEYS.LABELS_KEY].ravel().long()
            z_mean = inference_outputs[MODULE_KEYS.QZ_KEY].loc
            y_logits = self.classifier_0(z_mean)
            classification_loss += torch.nn.functional.cross_entropy(y_logits, y, reduction="mean")
        if self.classifier_1 is not None and mode == self.input_names[1]:
            y = tensors[REGISTRY_KEYS.LABELS_KEY].ravel().long()
            z_mean = inference_outputs[MODULE_KEYS.QZ_KEY].loc
            y_logits = self.classifier_1(z_mean)
            classification_loss += torch.nn.functional.cross_entropy(y_logits, y, reduction="mean")

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local_norm,
            extra_metrics={
                "z": inference_outputs[MODULE_KEYS.Z_KEY],
                "mu_all": mu_all,
                "logvar_all": logvar_all,
                "v_all": v_all,
                "guidance_graph": self.guidance_graph,
                "classification_loss": classification_loss,
            },
        )
