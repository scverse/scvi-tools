from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import torch
from torch.distributions import (
    Categorical,
    Exponential,
    Independent,
    MixtureSameFamily,
    Normal,
)
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, EmbeddingModuleMixin, LossOutput, auto_move_data
from scvi.nn import Encoder, FCLayers

if TYPE_CHECKING:
    from collections import OrderedDict
    from typing import Literal

    import numpy as np


def identity(x):
    """Identity function."""
    return x


class MRDeconv(EmbeddingModuleMixin, BaseModuleClass):
    """Model for multi-resolution deconvolution of spatial transriptomics.

    Parameters
    ----------
    n_spots
        Number of input spots
    n_labels
        Number of cell types
    n_hidden
        Number of neurons in the hidden layers
    n_layers
        Number of layers used in the encoder networks
    n_latent
        Number of dimensions used in the latent variables
    n_genes
        Number of genes used in the decoder
    px_r
        parameters for the px_r tensor in the CondSCVI model
    per_ct_bias
        estimates of per cell-type expression bias in the CondSCVI model
    decoder_state_dict
        state_dict from the decoder of the CondSCVI model
    px_decoder_state_dict
        state_dict from the px_decoder of the CondSCVI model
    dropout_decoder
        Dropout rate for the decoder neural network (same dropout as in CondSCVI decoder)
    dropout_amortization
        Dropout rate for the amortization neural network
    n_samples_augmentation
        Number of samples used in the augmentation
    n_states_per_label
        Number of states per cell-type in each spot
    eps_v
        Epsilon value for each cell-type proportion used during training.
    n_states_per_augmented_label
        Number of states per cell-type in each spot during augmentation
    mean_vprior
        Mean parameter for each component in the empirical prior over the latent space
    var_vprior
        Diagonal variance parameter for each component in the empirical prior over the latent space
    mp_vprior
        Mixture proportion in cell type sub-clustering of each component in the empirical prior
        amortization
    prior_mode
        Mode of the prior distribution for the latent space.
        Either "mog" for mixture of gaussians or "normal" for normal distribution.
    add_celltypes
        Number of additional cell types compared to single cell data to add to the model
    n_latent_amortization
        Number of dimensions used in the latent variables for the amortization encoder
        neural network
    extra_encoder_kwargs
        Extra keyword arguments passed into :class:`~scvi.nn.FCLayers`.
    extra_decoder_kwargs
        Extra keyword arguments passed into :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_spots: int,
        n_labels: int,
        n_batch: int,
        n_hidden: int,
        n_layers: int,
        n_latent: int,
        n_genes: int,
        decoder_state_dict: OrderedDict,
        px_decoder_state_dict: OrderedDict,
        px_r: torch.tensor,
        per_ct_bias: torch.tensor,
        dropout_decoder: float,
        dropout_amortization: float = 0.03,
        augmentation: bool = True,
        n_samples_augmentation: int = 2,
        n_states_per_label: int = 3,
        eps_v: float = 2e-3,
        mean_vprior: np.ndarray = None,
        var_vprior: np.ndarray = None,
        mp_vprior: np.ndarray = None,
        amortization: Literal["none", "latent", "proportion", "both"] = "both",
        prior_mode: Literal["mog", "normal"] = "mog",
        add_celltypes: int = 1,
        n_latent_amortization: int | None = None,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
    ):
        super().__init__()
        if prior_mode == "mog":
            assert amortization in ["both", "latent"], (
                "Amortization must be active for latent variables to use mixture "
                "of gaussians generation"
            )
        self.n_spots = n_spots
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.augmentation = augmentation
        self.n_samples_augmentation = n_samples_augmentation
        self.dropout_decoder = dropout_decoder
        self.n_states_per_label = n_states_per_label
        self.dropout_amortization = dropout_amortization
        self.n_genes = n_genes
        self.amortization = amortization
        self.prior_mode = prior_mode
        self.add_celltypes = add_celltypes
        self.eps_v = eps_v
        self.n_latent_amortization = n_latent_amortization
        # unpack and copy parameters
        _extra_decoder_kwargs = extra_decoder_kwargs or {}
        cat_list = [n_labels]
        self.init_embedding(REGISTRY_KEYS.BATCH_KEY, n_batch)
        batch_dim = self.get_embedding(REGISTRY_KEYS.BATCH_KEY).embedding_dim
        n_input_decoder = n_latent + batch_dim

        self.decoder = FCLayers(
            n_in=n_input_decoder,
            n_out=n_hidden,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_decoder,
            use_layer_norm=True,
            use_batch_norm=False,
            **_extra_decoder_kwargs,
        )
        self.decoder.load_state_dict(decoder_state_dict)
        for param in self.decoder.parameters():
            param.requires_grad = False
        self.px_decoder = torch.nn.Linear(n_hidden, n_genes)
        self.px_decoder.load_state_dict(px_decoder_state_dict)
        for param in self.px_decoder.parameters():
            param.requires_grad = False
        self.px_r = torch.nn.Parameter(px_r)
        self.register_buffer("per_ct_bias", per_ct_bias)

        # cell_type specific factor loadings
        self.V = torch.nn.Parameter(torch.randn(self.n_labels + self.add_celltypes, self.n_spots))

        # within cell_type factor loadings
        self.gamma = torch.nn.Parameter(torch.randn(n_latent, self.n_labels, self.n_spots))
        if mean_vprior is not None:
            self.register_buffer("mean_vprior", mean_vprior)
            self.register_buffer("var_vprior", var_vprior)
            self.register_buffer("mp_vprior", mp_vprior)
            cats = Categorical(probs=self.mp_vprior)
            normal_dists = Independent(
                Normal(self.mean_vprior, torch.sqrt(self.var_vprior) + 1e-4),
                reinterpreted_batch_ndims=1,
            )
            self.qz_prior = MixtureSameFamily(cats, normal_dists)
        else:
            self.mean_vprior = None
            self.var_vprior = None
        # noise from data
        self.eta = torch.nn.Parameter(torch.zeros(self.add_celltypes, self.n_genes))
        # additive gene bias
        self.beta = torch.nn.Parameter(torch.zeros(self.n_genes))

        # create additional neural nets for amortization
        # within cell_type factor loadings
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        if self.prior_mode == "mog":
            return_dist = (
                self.n_states_per_label * n_labels * n_latent + self.n_states_per_label * n_labels
            )
        else:
            return_dist = n_labels * n_latent
        if self.n_latent_amortization is not None:
            # Uses a combined latent space for proportions and gammas.
            self.z_encoder = Encoder(
                self.n_genes,
                n_latent_amortization,
                n_cat_list=[n_batch],
                n_layers=n_layers,
                n_hidden=n_hidden,
                dropout_rate=dropout_amortization,
                inject_covariates=True,
                use_batch_norm=False,
                use_layer_norm=True,
                var_activation=torch.nn.functional.softplus,
                return_dist=True,
                **_extra_encoder_kwargs,
            )
        else:

            def identity(x, batch_index=None):
                return x, Normal(x, scale=1e-6 * torch.ones_like(x))

            self.z_encoder = identity
            n_latent_amortization = self.n_genes
        self.gamma_encoder = torch.nn.Sequential(
            FCLayers(
                n_in=n_latent_amortization,
                n_out=n_hidden,
                n_cat_list=None,
                n_layers=2,
                n_hidden=n_hidden,
                dropout_rate=dropout_amortization,
                use_layer_norm=True,
                use_batch_norm=False,
            ),
            torch.nn.Linear(n_hidden, return_dist),
        )
        # cell type loadings
        self.V_encoder = torch.nn.Sequential(
            FCLayers(
                n_in=n_latent_amortization,
                n_out=n_hidden,
                n_cat_list=None,
                n_layers=n_layers,
                n_hidden=n_hidden,
                dropout_rate=0,
                use_layer_norm=True,
                use_batch_norm=False,
            ),
            torch.nn.Linear(n_hidden, n_labels + self.add_celltypes),
        )

    def _get_inference_input(self, tensors):
        x = tensors[REGISTRY_KEYS.X_KEY]
        x_smoothed = tensors.get("x_smoothed", None)
        m = x.shape[0]
        if x_smoothed is not None:
            n_samples = self.n_samples_augmentation + 2
            n_samples_observed = 2
        else:
            n_samples = self.n_samples_augmentation + 1
            n_samples_observed = 1
        px_r = torch.exp(self.px_r)
        if self.augmentation and self.training:
            with torch.no_grad():
                prior_sampled = self.qz_prior.sample([n_samples, m]).reshape(
                    n_samples, -1, self.n_latent
                )
                enum_label = (
                    torch.arange(0, self.n_labels).repeat(m).view((-1, 1))
                )  # m * n_labels, 1
                batch_rep = self.compute_embedding(
                    REGISTRY_KEYS.BATCH_KEY, tensors[REGISTRY_KEYS.BATCH_KEY]
                )
                batch_rep_input = batch_rep.repeat_interleave(self.n_labels, dim=0).repeat(
                    n_samples, 1, 1
                )
                decoder_input = torch.cat([prior_sampled, batch_rep_input], dim=-1)
                px_scale_augment_ = torch.nn.Softmax(dim=-1)(
                    self.px_decoder(self.decoder(decoder_input, enum_label.to(x.device)))
                    + self.per_ct_bias[enum_label.ravel()].unsqueeze(-3)
                    + self.beta.view(1, 1, -1)
                )
                px_scale_augment = px_scale_augment_.reshape(
                    (n_samples, x.shape[0], self.n_labels, -1)
                )  # (samples, mi, n_labels, n_genes)
                library = x.sum(-1).view(1, m, 1, 1).repeat(n_samples, 1, 1, 1)
                library[1, ...] = library[1, ...] + 50
                px_scale_augment = px_scale_augment.reshape(
                    n_samples, m, self.n_labels, -1
                )  # (samples, m, n_labels, n_genes)
                px_rate = library * px_scale_augment  # (samples, m, n_labels, n_genes)
                ratios_ct_augmentation = (
                    torch.distributions.Dirichlet(torch.zeros(self.n_labels) + 0.03)
                    .sample([n_samples, m])
                    .to(x.device)
                )
                ratios_ct_augmentation = ratios_ct_augmentation.reshape(
                    n_samples, m, self.n_labels
                )
                # sum over celltypes
                augmentation_rate = torch.einsum(
                    "imc, imcg -> img", ratios_ct_augmentation, px_rate
                )  # (samples, m, n_genes)
                ratio_augmentation_ = (
                    torch.distributions.Beta(0.4, 0.5)
                    .sample([self.n_samples_augmentation - 1, m])
                    .unsqueeze(-1)
                    .to(x.device)
                )
                ratio_augmentation = torch.cat(
                    [
                        torch.zeros((n_samples_observed, m, 1), device=x.device),
                        torch.ones((1, m, 1), device=x.device),
                        ratio_augmentation_,
                    ],
                    dim=0,
                )
                augmented_counts = NegativeBinomial(
                    mu=augmentation_rate, theta=px_r
                ).sample()  # (samples, m, n_genes)
                x_augmented = (1 - ratio_augmentation) * x.unsqueeze(
                    0
                ) + ratio_augmentation * augmented_counts
                if x_smoothed is not None:
                    x_augmented[1, ...] = x_smoothed
        else:
            x_augmented = x.unsqueeze(0)
            if x_smoothed is not None:
                x_augmented = torch.cat([x_augmented, x_smoothed.unsqueeze(0)], dim=0)
            prior_sampled = None
            ratios_ct_augmentation = None
            ratio_augmentation = None

        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            "batch_index": batch_index,
            "x_augmented": x_augmented,
            "prior_sampled": prior_sampled,
            "ratios_ct_augmentation": ratios_ct_augmentation,
            "ratio_augmentation": ratio_augmentation,
        }
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        ind_x = tensors[REGISTRY_KEYS.INDICES_KEY].long().ravel()
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {"z": z, "ind_x": ind_x, "library": library, "batch_index": batch_index}
        return input_dict

    @auto_move_data
    def inference(
        self,
        x_augmented,
        batch_index,
        n_samples=1,
        prior_sampled=None,
        ratios_ct_augmentation=None,
        ratio_augmentation=None,
    ):
        """Runs the inference (encoder) model."""
        x_ = x_augmented
        library = x_augmented.sum(-1).unsqueeze(-1)
        x_ = torch.log(1 + x_)
        if self.n_latent_amortization is not None:
            qz, z = self.z_encoder(x_, batch_index)
        else:
            z = x_
            qz = Normal(x_, scale=1e-6 * torch.ones_like(x_))

        outputs = {
            "z": z,
            "qz": qz,
            "library": library,
            "x_augmented": x_augmented,
            "prior_sampled": prior_sampled,
            "ratio_augmentation": ratio_augmentation,
            "ratios_ct_augmentation": ratios_ct_augmentation,
        }
        return outputs

    @auto_move_data
    def generative(self, z, ind_x, library, batch_index):
        """Build the deconvolution model for every cell in the minibatch."""
        m = len(ind_x)
        # setup all non-linearities
        eps = torch.nn.functional.softmax(self.eta, dim=-1)  # add_celltypes, n_genes
        px_r = torch.exp(self.px_r)
        n_samples = z.size(0)

        if self.amortization in ["both", "latent"]:
            if self.prior_mode == "mog":
                gamma_ = self.gamma_encoder(z)
                proportion_modes_logits = (
                    torch.transpose(
                        gamma_[:, :, -self.n_states_per_label * self.n_labels :], -1, -2
                    )
                    .reshape((n_samples, self.n_states_per_label, self.n_labels, m))
                    .transpose(-1, -2)
                )  # n_samples, n_states_per_label, m, n_labels
                proportion_modes = torch.nn.functional.softmax(proportion_modes_logits, dim=-3)
                gamma_ind = torch.transpose(
                    gamma_[:, :, : self.n_states_per_label * self.n_labels * self.n_latent], -1, -2
                ).reshape(
                    (n_samples, self.n_states_per_label, self.n_latent, self.n_labels, m)
                )  # n_samples, n_states_per_label, n_latent, n_labels, m
            else:
                gamma_ind = torch.transpose(self.gamma_encoder(z), 0, 1).reshape(
                    (n_samples, 1, self.n_latent, self.n_labels, m)
                )
                proportion_modes_logits = proportion_modes = torch.ones(
                    (n_samples, 1, m, self.n_labels), device=z.device
                )
        else:
            gamma_ind = (
                self.gamma[:, :, ind_x].unsqueeze(0).unsqueeze(0).repeat(n_samples, 1, 1, 1, 1)
            )  # n_samples, n_latent, n_labels, m
            proportion_modes_logits = proportion_modes = torch.ones(
                (n_samples, 1, m, self.n_labels), device=z.device
            )

        if self.amortization in ["both", "proportion"]:
            v_ind = self.V_encoder(z)
        else:
            v_ind = (
                self.V[:, ind_x].T.unsqueeze(0).repeat(n_samples, 1, 1)
            )  # n_samples, m, labels + 1
        v_ind = torch.nn.functional.softmax(v_ind, dim=-1)

        px_est = torch.zeros((n_samples, m, self.n_labels, self.n_genes), device=z.device)
        enum_label = torch.arange(0, self.n_labels).repeat(m).view((-1, 1))  # m * n_labels, 1

        batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
        batch_rep_input = (
            batch_rep.repeat_interleave(self.n_labels, dim=0).unsqueeze(0).repeat(n_samples, 1, 1)
        )

        for mode in range(gamma_ind.shape[-4]):
            # reshape and get gene expression value for all minibatch
            gamma_ind_ = torch.transpose(
                gamma_ind[:, mode, ...], -1, -3
            )  # n_samples, m, n_labels, n_latent
            gamma_reshape_ = gamma_ind_.reshape(
                (n_samples, -1, self.n_latent)
            )  # n_samples, m * n_labels, n_latent
            decoder_input_ = torch.cat([gamma_reshape_, batch_rep_input], dim=-1)
            h = self.decoder(decoder_input_, enum_label.to(z.device))
            px_est += proportion_modes[:, mode, ...].unsqueeze(-1) * torch.nn.Softmax(dim=-1)(
                self.px_decoder(h).reshape((n_samples, m, self.n_labels, -1))
                + self.per_ct_bias[enum_label.ravel()].reshape(1, m, self.n_labels, -1)
                + self.beta.view(1, 1, 1, -1)
            )  # (n_samples, m, n_labels, n_genes)

        # add the additional cell types
        eps = eps.unsqueeze(0).repeat(
            n_samples, m, 1, 1
        )  # (n_samples, m, add_celltypes, n_genes) <- this is the dummy cell type
        r_hat = torch.cat([px_est, eps], dim=-2)  # n_samples, m, n_labels + add_celltypes, n_genes

        # now combine them for convolution. Add epsilon during training.
        eps_v = self.eps_v if self.training else 0.0
        px_scale = torch.sum(
            (v_ind.unsqueeze(-1) + eps_v) * r_hat, dim=-2
        )  # n_samples, m, n_genes
        px_rate = library * px_scale

        return {
            "px_r": px_r,
            "px_rate": px_rate,
            "px_scale": px_scale,
            "px_mu": r_hat,
            "gamma": gamma_ind,
            "v": v_ind,
            "proportion_modes": proportion_modes,
            "proportion_modes_logits": proportion_modes_logits,
        }

    def _compute_cross_entropy(self, prob_true, prob_pred):
        log_prob_pred = torch.log(prob_pred / prob_pred.sum(axis=-1, keepdim=True))
        prob_true = prob_true + 1e-20
        prob_true = prob_true / prob_true.sum(axis=-1, keepdim=True)
        kl_div = torch.nn.functional.kl_div(
            log_prob_pred, prob_true, reduction="batchmean", log_target=False
        )

        return kl_div

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        ct_sparsity_weight: float = 2.0,
        weighting_augmentation: float = 100.0,
        weighting_smoothing: float = 100.0,
        eta_reg: float = 1.0,
        beta_reg: float = 1.0,
        weighting_kl_latent: float = 1e-3,
        reconst_weight: float = 3.0,
    ):
        """Compute the loss."""
        x_augmented = inference_outputs["x_augmented"]
        px_rate = generative_outputs["px_rate"]
        px_r = generative_outputs["px_r"]
        gamma = generative_outputs["gamma"]
        v = generative_outputs["v"]
        ratio_augmentation = inference_outputs["ratio_augmentation"]
        ratios_ct_augmentation = inference_outputs["ratios_ct_augmentation"]

        if "x_smoothed" in tensors:
            sample_fully_augmented = 2
            n_samples = self.n_samples_augmentation + 2
        else:
            sample_fully_augmented = 1
            n_samples = self.n_samples_augmentation + 1
        m = x_augmented.shape[1]

        if self.augmentation:
            prior_sampled = inference_outputs["prior_sampled"].reshape(
                n_samples, 1, x_augmented.shape[1], self.n_labels, self.n_latent
            )
            mean_vprior = torch.cat(
                [
                    self.mean_vprior.unsqueeze(0).unsqueeze(-4).repeat(n_samples, m, 1, 1, 1),
                    prior_sampled.permute(0, 2, 3, 1, 4),
                ],
                dim=-2,
            )  # n_samples, m, n_labels, p + 1, n_latent
            var_vprior = torch.cat(
                [
                    self.var_vprior.unsqueeze(0).unsqueeze(-4).repeat(n_samples, m, 1, 1, 1),
                    torch.min(self.var_vprior, dim=-2)
                    .values.view(1, 1, self.n_labels, 1, self.n_latent)
                    .repeat(n_samples, m, 1, 1, 1),
                ],
                dim=-2,
            )  # n_samples, m, n_labels, p + 1, n_latent

            mp_vprior = torch.cat(
                [
                    (1 - ratio_augmentation.unsqueeze(-1))
                    * self.mp_vprior.view(1, 1, self.n_labels, -1).repeat(n_samples, m, 1, 1),
                    ratio_augmentation.unsqueeze(-1) * ratios_ct_augmentation.unsqueeze(-1),
                ],
                dim=-1,
            )  # n_samples, m, n_labels, p + 1
        else:
            mean_vprior = self.mean_vprior.unsqueeze(0).unsqueeze(0)
            var_vprior = self.var_vprior.unsqueeze(0).unsqueeze(0)
            mp_vprior = self.mp_vprior.unsqueeze(0).unsqueeze(0)

        proportion_modes = generative_outputs["proportion_modes"]
        # rounding errors for softmax lead to px_rate 0 which induces NaNs in the log_prob
        reconst_loss = (
            -NegativeBinomial(mu=px_rate, theta=px_r).log_prob(x_augmented).sum(-1).mean(0)
        )
        # beta prior likelihood
        mean = torch.zeros_like(self.beta)
        scale = torch.ones_like(self.beta)
        # beta loss
        glo_neg_log_likelihood_prior = -beta_reg * Normal(mean, scale).log_prob(self.beta).sum()
        loss_augmentation = torch.tensor(0.0, device=x_augmented.device)
        if self.augmentation:
            expected_proportions = torch.cat(
                [
                    ratios_ct_augmentation,
                    torch.zeros([n_samples, m, self.add_celltypes]).to(v.device),
                ],
                dim=-1,
            )
            loss_augmentation += weighting_augmentation * self._compute_cross_entropy(
                expected_proportions[sample_fully_augmented, ...].squeeze(0),
                v[sample_fully_augmented, ...].squeeze(0),
            )
        if "x_smoothed" in tensors:
            loss_augmentation += weighting_smoothing * self._compute_cross_entropy(
                v[1, ...].squeeze(0), v[0, ...].squeeze(0)
            )

        # gamma prior likelihood
        if self.mean_vprior is None:
            # isotropic normal prior
            mean = torch.zeros_like(gamma)
            scale = torch.ones_like(gamma)
            neg_log_likelihood_prior = -Normal(mean, scale).log_prob(gamma).sum(2).sum(1)
        elif self.prior_mode == "mog":
            # gamma is of shape n_samples, minibatch_size, 1, n_latent, n_labels
            gamma = gamma.permute(
                1, 0, 4, 3, 2
            )  # p, n_samples, minibatch_size, n_labels, n_latent
            cats = Categorical(probs=mp_vprior)
            normal_dists = Independent(
                Normal(mean_vprior, var_vprior), reinterpreted_batch_ndims=1
            )
            pre_lse = MixtureSameFamily(cats, normal_dists).log_prob(
                gamma
            )  # p, n_samples, minibatch_size, n_labels
            pre_lse = pre_lse.permute(1, 0, 2, 3)
            log_likelihood_prior = torch.mul(pre_lse, proportion_modes + 1e-6).sum(
                -3
            )  # n_samples, minibatch, n_labels
            neg_log_likelihood_prior = -torch.mul(
                log_likelihood_prior, v[:, :, : -self.add_celltypes] + 1e-3
            ).sum(-1)  # n_samples, minibatch, n_labels
            # neg_log_likelihood_prior = - log_likelihood_prior.sum(-1)
        else:
            gamma = gamma.permute(
                0, 4, 1, 3, 2
            )  # n_samples, minibatch_size, 1, n_labels, n_latent
            mean_vprior = torch.transpose(
                mean_vprior, -3, -2
            )  # n_samples, m, p, n_labels, n_latent
            var_vprior = torch.transpose(var_vprior, -3, -2)  # n_samples, m, p, n_labels, n_latent
            mp_vprior = torch.transpose(mp_vprior, -2, -1)  # n_samples, m, p, n_labels
            pre_lse = (
                Normal(mean_vprior, torch.sqrt(var_vprior) + 1e-4)
                .log_prob(gamma)
                .sum(-1)
                .squeeze(-3)
            ) + torch.log(1e-3 + mp_vprior)  # n_samples, minibatch, p, n_labels
            log_likelihood_prior = torch.logsumexp(pre_lse, -2)  # n_samples, minibatch, n_labels
            neg_log_likelihood_prior = -log_likelihood_prior.sum(-1)  # n_samples, minibatch
        if self.n_latent_amortization is not None:
            neg_log_likelihood_prior += weighting_kl_latent * kl(
                inference_outputs["qz"],
                Normal(
                    torch.zeros([self.n_latent_amortization], device=x_augmented.device),
                    torch.ones([self.n_latent_amortization], device=x_augmented.device),
                ),
            ).sum(dim=-1)

        v_sparsity_loss = (
            ct_sparsity_weight * torch.distributions.Categorical(probs=v[0, ...]).entropy()
        )
        v_sparsity_loss -= (
            eta_reg
            * Exponential(torch.ones_like(v[0, :, -self.add_celltypes :]))
            .log_prob(v[0, :, -self.add_celltypes :])
            .sum()
        )

        loss = torch.mean(
            reconst_weight * reconst_loss
            + kl_weight * (neg_log_likelihood_prior + v_sparsity_loss)
            + glo_neg_log_likelihood_prior
            + loss_augmentation
        )

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=neg_log_likelihood_prior,
            kl_global=glo_neg_log_likelihood_prior,
            extra_metrics={
                "v_sparsity": v_sparsity_loss.mean(),
                "augmentation": loss_augmentation.mean(),
            },
        )

    @torch.inference_mode()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        """Sample from the posterior."""
        raise NotImplementedError("No sampling method for DestVI")
