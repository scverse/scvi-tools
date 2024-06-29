from collections import OrderedDict
from typing import Literal, Optional

import numpy as np
import torch
from torch.distributions import Categorical, Independent, MixtureSameFamily, Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi._types import Tunable
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder, FCLayers


def identity(x):
    """Identity function."""
    return x


class MRDeconv(BaseModuleClass):
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
    dropout_decoder
        Dropout rate for the decoder neural network (same dropout as in CondSCVI decoder)
    dropout_amortization
        Dropout rate for the amortization neural network
    decoder_state_dict
        state_dict from the decoder of the CondSCVI model
    px_decoder_state_dict
        state_dict from the px_decoder of the CondSCVI model
    px_r
        parameters for the px_r tensor in the CondSCVI model
    mean_vprior
        Mean parameter for each component in the empirical prior over the latent space
    var_vprior
        Diagonal variance parameter for each component in the empirical prior over the latent space
    mp_vprior
        Mixture proportion in cell type sub-clustering of each component in the empirical prior
        amortization
    beta_reg
        Scalar parameter indicating the strength of the variance penalty for
        the multiplicative offset in gene expression values (beta parameter). Default is 5
        (setting to 0.5 might help if single cell reference and spatial assay are different
        e.g. UMI vs non-UMI.)
    eta_reg
        Scalar parameter indicating the strength of the prior for
        the noise term (eta parameter). Default is 1e-4.
        (changing value is discouraged.)
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
        n_hidden: Tunable[int],
        n_layers: Tunable[int],
        n_latent: Tunable[int],
        n_genes: int,
        decoder_state_dict: OrderedDict,
        px_decoder_state_dict: OrderedDict,
        px_r: np.ndarray,
        dropout_decoder: float,
        augmentation: bool = False,
        n_samples_augmentation: int = 1,
        n_states_per_label: Tunable[int] = 1,
        n_states_per_augmented_label: float = 1,
        dropout_amortization: float = 0.03,
        mean_vprior: np.ndarray = None,
        var_vprior: np.ndarray = None,
        mp_vprior: np.ndarray = None,
        amortization: Literal["none", "latent", "proportion", "both"] = "both",
        beta_reg: Tunable[float] = 500.0,
        eta_reg: Tunable[float] = 1e-7,
        prior_mode: Literal["mog", "normal"] = "normal",
        n_latent_amortization: Optional[int] = None,
        extra_encoder_kwargs: Optional[dict] = None,
        extra_decoder_kwargs: Optional[dict] = None,
    ):
        super().__init__()
        if prior_mode == 'mog':
            assert amortization in ["both", "latent"], "Amortization must be active for latent variables to use mixture of gaussians generation"
        self.n_spots = n_spots
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.augmentation = augmentation
        self.n_samples_augmentation = n_samples_augmentation
        self.n_states_per_augmented_label = n_states_per_augmented_label
        self.dropout_decoder = dropout_decoder
        self.n_states_per_label = n_states_per_label
        self.dropout_amortization = dropout_amortization
        self.n_genes = n_genes
        self.amortization = amortization
        self.beta_reg = beta_reg
        self.eta_reg = eta_reg
        self.prior_mode = prior_mode
        self.n_latent_amortization = n_latent_amortization
        # unpack and copy parameters
        _extra_decoder_kwargs = extra_decoder_kwargs or {}
        cat_list = [n_labels]
        self.init_embedding(REGISTRY_KEYS.BATCH_KEY, n_batch, {})
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
        self.px_o = torch.nn.Parameter(px_r)

        # cell_type specific factor loadings
        self.V = torch.nn.Parameter(torch.randn(self.n_labels + 1, self.n_spots))

        # within cell_type factor loadings
        self.gamma = torch.nn.Parameter(
            torch.randn(n_latent, self.n_labels, self.n_spots)
        )
        if mean_vprior is not None:
            self.register_buffer("mean_vprior", mean_vprior)
            self.register_buffer("var_vprior", var_vprior)
            self.register_buffer("mp_vprior", mp_vprior)
            cats = Categorical(probs=self.mp_vprior)
            normal_dists = Independent(
                Normal(
                    self.mean_vprior,
                    torch.sqrt(self.var_vprior) + 1e-4
                ),
                reinterpreted_batch_ndims=1
            )
            self.qz_prior = MixtureSameFamily(cats, normal_dists)
        else:
            self.mean_vprior = None
            self.var_vprior = None
        # noise from data
        self.eta = torch.nn.Parameter(torch.zeros(self.n_genes))
        # additive gene bias
        self.beta = torch.nn.Parameter(torch.zeros(self.n_genes))
        print('beta is parameter')

        # create additional neural nets for amortization
        # within cell_type factor loadings
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        if self.prior_mode == "mog":
            return_dist = self.n_states_per_label * n_labels * n_latent + self.n_states_per_label * n_labels
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
                return x, Normal(x, scale=1e-6*torch.ones_like(x))
            self.z_encoder = identity
            n_latent_amortization = self.n_genes
            n_layers = 2
        self.gamma_encoder = torch.nn.Sequential(
            FCLayers(
                n_in=n_latent_amortization,
                n_out=n_hidden,
                n_cat_list=None,
                n_layers=n_layers,
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
            torch.nn.Linear(n_hidden, n_labels + 1),
        )

    def _get_inference_input(self, tensors):
        x = tensors[REGISTRY_KEYS.X_KEY]
        m = x.shape[0]
        n_samples = self.n_samples_augmentation + 1
        if self.augmentation and self.training:
            with torch.no_grad():
                beta = torch.exp(self.beta)  # n_genes
                # beta = torch.cat([beta.view(1, 1, 1, -1), torch.ones_like(beta).view(1, 1, 1, -1).repeat(n_samples-1, 1, 1, 1)])
                prior_sampled = self.qz_prior.sample(
                    [n_samples, self.n_states_per_augmented_label, m]).reshape(
                        n_samples*self.n_states_per_augmented_label, -1, self.n_latent)
                enum_label = (
                    torch.arange(0, self.n_labels).repeat(m).view((-1, 1))
                )  # m * n_labels, 1
                batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, tensors[REGISTRY_KEYS.BATCH_KEY])
                batch_rep_input = batch_rep.repeat_interleave(self.n_labels, dim=0)
                decoder_input = torch.cat([prior_sampled, batch_rep_input], dim=-1)
                px_scale_augment_ = torch.nn.Softmax(dim=-1)(self.px_decoder(self.decoder(decoder_input, enum_label.to(x.device))) + beta.view(1, 1, -1))
                px_scale_augment = px_scale_augment_.reshape(
                    (n_samples*self.n_states_per_augmented_label, x.shape[0], self.n_labels, -1)
                )  # (samples*states_per_cell, mi, n_labels, n_genes)
                library = x.sum(-1).view(1, 1, m, 1, 1).repeat(n_samples, 1, 1, 1, 1)
                library[1, ...] = library[1, ...] + 50
                px_scale_augment = px_scale_augment.reshape(n_samples, self.n_states_per_augmented_label, m, self.n_labels, -1) # (samples, states_per_cell, m, n_labels, n_genes)
                px_rate = library * px_scale_augment # (samples, states_per_cell, m, n_labels, n_genes)
                ratios_ct_augmentation = torch.distributions.Dirichlet(
                    torch.zeros(self.n_states_per_augmented_label * self.n_labels) + 0.03).sample([n_samples, m]).to(x.device)
                ratios_ct_augmentation = ratios_ct_augmentation.reshape(n_samples, m, self.n_states_per_augmented_label, self.n_labels).permute(0, 2, 1, 3)
                augmentation_rate = torch.einsum('ilmk, ilmkg -> img', ratios_ct_augmentation, px_rate)  # (samples, m, n_genes)
                ratio_augmentation_ = torch.distributions.Beta(0.4, 0.5).sample([self.n_samples_augmentation-1, m]).unsqueeze(-1).to(x.device)
                ratio_augmentation = torch.cat([torch.zeros((1, m, 1), device=x.device), torch.ones((1, m, 1), device=x.device), ratio_augmentation_], dim=0)
                augmented_counts = NegativeBinomial(
                    augmentation_rate, logits=self.px_o
                ).sample()  # (samples*states_per_cell, m, n_labels, n_genes)
                # print('TTTT1', augmentation_rate[1, ...].sum(-1).min(), augmented_counts[1, ...].sum(-1).min(), x.sum(-1).min(), x.shape)
                x_augmented = (
                    (1 - ratio_augmentation) * x.unsqueeze(0) +
                    ratio_augmentation * augmented_counts
                )
        else:
            x_augmented = x.unsqueeze(0)
            prior_sampled = None
            ratios_ct_augmentation = None
            ratio_augmentation = None
            
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            "batch_index": batch_index,
            "x_augmented": x_augmented,
            "prior_sampled": prior_sampled,
            "ratios_ct_augmentation": ratios_ct_augmentation,
            "ratio_augmentation": ratio_augmentation}
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
            qz = Normal(x_, scale=1e-6*torch.ones_like(x_))

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
        beta = torch.exp(self.beta)  # n_genes
        eps = torch.nn.functional.softmax(self.eta, dim=-1) # n_genes
        if self.training and self.augmentation:
            n_samples = (self.n_samples_augmentation + 1)
        else:
            n_samples = 1

        if self.amortization in ["both", "latent"]:
            if self.prior_mode == "mog":
                gamma_ = self.gamma_encoder(z)
                proportion_modes_logits = torch.transpose(
                    gamma_[:, :, -self.n_states_per_label*self.n_labels:], -1, -2).reshape(
                    (n_samples, self.n_states_per_label, self.n_labels, m)
                ).transpose(-1, -2) # n_samples, n_states_per_label, m, n_labels
                proportion_modes = torch.nn.functional.softmax(proportion_modes_logits, dim=-3)
                gamma_ind = torch.transpose(
                    gamma_[:, :, :self.n_states_per_label*self.n_labels*self.n_latent], -1, -2).reshape(
                        (n_samples, self.n_states_per_label, self.n_latent, self.n_labels, m)
                ) # n_samples, n_states_per_label, n_latent, n_labels, m
            else:
                gamma_ind = torch.transpose(
                    self.gamma_encoder(z), 0, 1).reshape(
                        (n_samples, 1, self.n_latent, self.n_labels, m)
                )
                proportion_modes_logits = proportion_modes = torch.ones(
                    (n_samples, 1, m, self.n_labels), device=z.device)
        else:
            gamma_ind = self.gamma[:, :, ind_x].unsqueeze(0).unsqueeze(0).repeat(
                n_samples, 1, 1, 1, 1)  # n_samples, n_latent, n_labels, m
            proportion_modes_logits = proportion_modes = torch.ones((n_samples, 1, m, self.n_labels), device=z.device)

        if self.amortization in ["both", "proportion"]:
            v_ind = self.V_encoder(z)
        else:
            v_ind = self.V[:, ind_x].T.unsqueeze(0).repeat(
                n_samples, 1, 1) # n_samples, m, labels + 1
        v_ind = torch.nn.functional.softmax(v_ind, dim=-1)

        px_est = torch.zeros((n_samples, m, self.n_labels, self.n_genes), device=z.device)
        enum_label = (
            torch.arange(0, self.n_labels).repeat(m).view((-1, 1))
        )  # m * n_labels, 1
        
        batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
        batch_rep_input = batch_rep.repeat_interleave(self.n_labels, dim=0).unsqueeze(0).repeat(n_samples, 1, 1)

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
            px_est += proportion_modes[:, mode, ...].unsqueeze(-1) * torch.nn.Softmax(dim=-1)(self.px_decoder(h).reshape(
                (n_samples, m, self.n_labels, -1)
            )  + beta.view(1, 1, 1, -1)) # (n_samples, m, n_labels, n_genes)

        # add the dummy cell type
        eps = eps.unsqueeze(0).repeat(n_samples, m, 1).unsqueeze(-2) # (n_samples, m, 1, n_genes) <- this is the dummy cell type

        # account for gene specific bias and add noise, take sample without augmentation.
        r_hat = torch.cat(
            [beta.view(1, 1, 1, -1) * px_est, eps], dim=-2
        )  # n_samples, m, n_labels + 1, n_genes
        
        # now combine them for convolution
        px_scale = torch.sum(v_ind.unsqueeze(-1) * r_hat, dim=-2)  # n_samples, m, n_genes
        px_rate = library * px_scale
        px_mu = torch.exp(self.px_o) * r_hat

        return {
            "px_o": self.px_o,
            "px_rate": px_rate,
            "px_mu": px_mu,
            "px_scale": px_scale,
            "gamma": gamma_ind,
            "v": v_ind,
            "proportion_modes": proportion_modes,
            "proportion_modes_logits": proportion_modes_logits,
        }

    def _compute_cross_entropy(self, prob_true, prob_pred):
        log_prob_pred = torch.log(prob_pred / prob_pred.sum(axis=-1, keepdim=True))
        prob_true = prob_true + 1e-20
        prob_true = prob_true / prob_true.sum(axis=-1, keepdim=True)
        kl_div = torch.nn.functional.kl_div(log_prob_pred, prob_true, reduction='batchmean', log_target=False)

        return kl_div

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        ct_sparsity_weight: float = 0.,
        weighting_augmentation: float = 10.,
    ):
        """Compute the loss."""
        x_augmented = inference_outputs["x_augmented"]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        gamma = generative_outputs["gamma"]
        v = generative_outputs["v"]
        ratio_augmentation = inference_outputs['ratio_augmentation']
        ratios_ct_augmentation = inference_outputs['ratios_ct_augmentation']
        
        n_samples = self.n_samples_augmentation + 1
        m = x_augmented.shape[1]

        if self.augmentation:
            prior_sampled = inference_outputs['prior_sampled'].reshape(n_samples, self.n_states_per_augmented_label, x_augmented.shape[1], self.n_labels, self.n_latent)
            mean_vprior = torch.cat(
                [self.mean_vprior.unsqueeze(0).unsqueeze(-4).repeat(n_samples, m, 1, 1, 1),
                 prior_sampled.permute(0, 2, 3, 1, 4)],
                dim=-2)  # n_samples, m, n_labels, p+n_states_per_augmented_label, n_latent
            var_vprior = torch.cat(
                [self.var_vprior.unsqueeze(0).unsqueeze(-4).repeat(n_samples, m, 1, 1, 1),
                 torch.min(self.var_vprior, dim=-2).values.view(
                     1, 1, self.n_labels, 1, self.n_latent).repeat(n_samples, m, 1, self.n_states_per_augmented_label, 1)],
                dim=-2
            ) # n_samples, m, n_labels, p+n_states_per_augmented_label, n_latent
            mp_vprior=torch.cat(
                [(1- ratio_augmentation.unsqueeze(-1)) * self.mp_vprior.view(1, 1, self.n_labels, -1).repeat(n_samples, m, 1, 1),
                 ratio_augmentation.unsqueeze(-1) * ratios_ct_augmentation.permute(0, 2, 3, 1)
                ],
                dim=-1
            ) # n_samples, m, n_labels, p+n_states_per_augmented_label
        else:
            mean_vprior = self.mean_vprior.unsqueeze(0).unsqueeze(0) 
            var_vprior = self.var_vprior.unsqueeze(0).unsqueeze(0) 
            mp_vprior = self.mp_vprior.unsqueeze(0).unsqueeze(0) 
 
        proportion_modes = generative_outputs["proportion_modes"]
        reconst_loss = - NegativeBinomial(px_rate, logits=px_o).log_prob(x_augmented).sum(-1)
        # eta prior likelihood
        mean = torch.zeros_like(self.eta)
        scale = torch.ones_like(self.eta)
        glo_neg_log_likelihood_prior = (
            -self.eta_reg * Normal(mean, scale).log_prob(self.eta).sum()
        )
        # beta loss
        glo_neg_log_likelihood_prior += (
            -self.beta_reg * Normal(mean, scale).log_prob(self.beta).sum()
        )
        if self.augmentation:
            expected_proportions = (
                ratio_augmentation * torch.cat([ratios_ct_augmentation.sum(-3), torch.zeros([n_samples, m, 1]).to(v.device)], dim=-1) +
                (1 - ratio_augmentation) * v[0, :, :] # unperturbed proportions
            )
            #loss_augmentation = weighting_augmentation * torch.abs(v - expected_proportions).sum(-1)
            loss_augmentation = 0
            for i in [1]:
                loss_augmentation += weighting_augmentation * self._compute_cross_entropy(expected_proportions[i, ...].squeeze(0), v[i, ...].squeeze(0))
        else:
            loss_augmentation = torch.tensor(0., device=x_augmented.device)

        # gamma prior likelihood
        if self.mean_vprior is None:
            # isotropic normal prior
            mean = torch.zeros_like(gamma)
            scale = torch.ones_like(gamma)
            neg_log_likelihood_prior = - Normal(mean, scale).log_prob(gamma).sum(2).sum(1)
        elif self.prior_mode == "mog":
            # gamma is of shape n_samples, minibatch_size, 1, n_latent, n_labels
            gamma = gamma.permute(1, 0, 4, 3, 2) # p, n_samples, minibatch_size, n_labels, n_latent
            cats = Categorical(probs=mp_vprior)
            normal_dists = Independent(
                Normal(
                    mean_vprior,
                    var_vprior
                ),
                reinterpreted_batch_ndims=1
            )
            pre_lse = MixtureSameFamily(cats, normal_dists).log_prob(gamma) # p, n_samples, minibatch_size, n_labels
            pre_lse = pre_lse.permute(1, 0, 2, 3)
            log_likelihood_prior = torch.mul(
                pre_lse,
                proportion_modes + 1e-3
            ).sum(-3) # n_samples, minibatch, n_labels
            neg_log_likelihood_prior = - log_likelihood_prior.sum(-1)
        else:
            gamma = gamma.permute(0, 4, 1, 3, 2) # n_samples, minibatch_size, 1, n_labels, n_latent
            mean_vprior = torch.transpose(mean_vprior, -3, -2) # n_samples, m, p, n_labels, n_latent
            var_vprior = torch.transpose(var_vprior, -3, -2) # n_samples, m, p, n_labels, n_latent
            mp_vprior = torch.transpose(mp_vprior, -2, -1)  # n_samples, m, p, n_labels
            pre_lse = (
                Normal(mean_vprior, torch.sqrt(var_vprior) + 1e-4).log_prob(gamma).sum(-1).squeeze(-3)
            ) + torch.log(1e-3 + mp_vprior)  # n_samples, minibatch, p, n_labels
            log_likelihood_prior = torch.logsumexp(pre_lse, -2)  # n_samples, minibatch, n_labels
            neg_log_likelihood_prior = - log_likelihood_prior.sum(-1)  # n_samples, minibatch
        if self.n_latent_amortization is not None:
            neg_log_likelihood_prior += kl(
                inference_outputs["qz"],
                Normal(torch.zeros([self.n_latent_amortization], device=x.device), torch.ones([self.n_latent_amortization], device=x.device))
            ).sum(dim=-1)
        
        v_sparsity_loss = ct_sparsity_weight * torch.distributions.Categorical(probs=v[0, :, :]).entropy()
        
        loss = torch.mean(
            reconst_loss + kl_weight * (neg_log_likelihood_prior + v_sparsity_loss) + glo_neg_log_likelihood_prior + loss_augmentation
        )

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=neg_log_likelihood_prior,
            kl_global=glo_neg_log_likelihood_prior,
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

