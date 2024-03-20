from collections import OrderedDict
from typing import Literal, Optional

import numpy as np
import torch
from torch.distributions import Normal
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
        which of the latent variables to amortize inference over (gamma, proportions, both or none)
    l1_reg
        Scalar parameter indicating the strength of L1 regularization on cell type proportions.
        A value of 50 leads to sparser results.
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
        n_batch_sc: Optional[list] = None,
        dropout_amortization: float = 0.05,
        mean_vprior: np.ndarray = None,
        var_vprior: np.ndarray = None,
        mp_vprior: np.ndarray = None,
        amortization: Literal["none", "latent", "proportion", "both"] = "both",
        l1_reg: Tunable[float] = 0.0,
        beta_reg: Tunable[float] = 5.0,
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
        self.dropout_decoder = dropout_decoder
        self.dropout_amortization = dropout_amortization
        self.n_genes = n_genes
        self.amortization = amortization
        self.l1_reg = l1_reg
        self.beta_reg = beta_reg
        self.eta_reg = eta_reg
        self.prior_mode = prior_mode
        self.n_latent_amortization = n_latent_amortization
        # unpack and copy parameters
        _extra_decoder_kwargs = extra_decoder_kwargs or {}
        cat_list = [n_labels, n_batch_sc]

        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_decoder,
            use_layer_norm=True,
            use_batch_norm=False,
            **_extra_decoder_kwargs,
        )
        self.px_decoder = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_genes), torch.nn.Softplus()
        )
        # don't compute gradient for those parameters
        self.decoder.load_state_dict(decoder_state_dict)
        for param in self.decoder.parameters():
            param.requires_grad = False
        self.px_decoder.load_state_dict(px_decoder_state_dict)
        for param in self.px_decoder.parameters():
            param.requires_grad = False
        self.register_buffer("px_o", torch.tensor(px_r))

        # cell_type specific factor loadings
        self.V = torch.nn.Parameter(torch.randn(self.n_labels + 1, self.n_spots))

        # within cell_type factor loadings
        self.gamma = torch.nn.Parameter(
            torch.randn(n_latent, self.n_labels, self.n_spots)
        )
        if mean_vprior is not None:
            self.p = mean_vprior.shape[1]
            self.register_buffer("mean_vprior", torch.tensor(mean_vprior))
            self.register_buffer("var_vprior", torch.tensor(var_vprior))
            self.register_buffer("mp_vprior", torch.tensor(mp_vprior))
        else:
            self.mean_vprior = None
            self.var_vprior = None
        # noise from data
        self.eta = torch.nn.Parameter(torch.randn(self.n_genes))
        # additive gene bias
        self.beta = torch.nn.Parameter(0.01 * torch.randn(self.n_genes))

        # create additional neural nets for amortization
        # within cell_type factor loadings
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        if self.prior_mode == "mog":
            print('Using mixture of gaussians for prior')
            return_dist = self.p * n_labels * n_latent + self.p * n_labels
        else:
            print('Using normal prior')
            return_dist = n_labels * n_latent
        print(f"return_dist: {return_dist}, {self.p}, {n_labels}, {n_latent}, {mean_vprior.shape}")
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
                dropout_rate=dropout_amortization,
                use_layer_norm=True,
                use_batch_norm=False,
            ),
            torch.nn.Linear(n_hidden, n_labels + 1),
        )

    def _get_inference_input(self, tensors):
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {"x": x, "batch_index": batch_index}
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        ind_x = tensors[REGISTRY_KEYS.INDICES_KEY].long().ravel()
        batch_index_sc = tensors['batch_index_sc'] + 2.

        input_dict = {"z": z, "ind_x": ind_x, "library": library, "batch_index_sc": batch_index_sc}
        return input_dict

    @auto_move_data
    def inference(
        self,
        x,
        batch_index,
        n_samples=1,
    ):
        """Runs the inference (encoder) model."""
        x_ = x
        library = torch.log(x.sum(1)).unsqueeze(1)
        x_ = torch.log(1 + x_)
        if self.n_latent_amortization is not None:
            qz, z = self.z_encoder(x_, batch_index)
        else:
            z = x_
            qz = Normal(x_, scale=1e-6*torch.ones_like(x_)) # dummy distribution

        outputs = {"z": z, "qz": qz, "library": library}
        return outputs

    @auto_move_data
    def generative(self, z, ind_x, library, batch_index_sc):
        """Build the deconvolution model for every cell in the minibatch."""
        m = len(ind_x)
        # setup all non-linearities
        beta = torch.exp(self.beta)  # n_genes
        eps = torch.nn.functional.softplus(self.eta) # n_genes

        if self.amortization in ["both", "latent"]:
            if self.prior_mode == "mog":
                gamma_ = self.gamma_encoder(z)
                proportion_modes_logits = torch.transpose(
                    gamma_[:, -self.p*self.n_labels:], 0, 1).reshape(
                    (self.p, self.n_labels, m)
                ).transpose(1, 2)
                proportion_modes = torch.nn.functional.softmax(proportion_modes_logits, dim=0)
                gamma_ind = torch.transpose(
                    gamma_[:, :self.p*self.n_labels*self.n_latent], 0, 1).reshape(
                        (self.p, self.n_latent, self.n_labels, -1)
                )
            else:
                gamma_ind = torch.transpose(
                    self.gamma_encoder(z), 0, 1).reshape(
                        (1, self.n_latent, self.n_labels, -1)
                )
                proportion_modes_logits = proportion_modes = torch.ones((1, self.n_labels), device=z.device)
        else:
            gamma_ind = self.gamma[:, :, ind_x].unsqueeze(0)  # 1, n_latent, n_labels, minibatch_size
            proportion_modes_logits = proportion_modes = torch.ones((1, self.n_labels), device=z.device)

        if self.amortization in ["both", "proportion"]:
            v_ind = self.V_encoder(z)
        else:
            v_ind = self.V[:, ind_x].T  # minibatch_size, labels + 1
        v_ind = torch.nn.functional.softplus(v_ind)

        px_est = torch.zeros((m, self.n_labels, self.n_genes), device=z.device)
        enum_label = (
            torch.arange(0, self.n_labels).repeat(m).view((-1, 1))
        )  # minibatch_size * n_labels, 1
        batch_index_sc_input = batch_index_sc.repeat_interleave(self.n_labels, dim=0)

        for mode in range(gamma_ind.shape[0]):
            # reshape and get gene expression value for all minibatch
            gamma_ind_ = torch.transpose(
                gamma_ind[mode, ...], 2, 0
            )  # minibatch_size, n_labels, n_latent
            gamma_reshape_ = gamma_ind_.reshape(
                (-1, self.n_latent)
            )  # minibatch_size * n_labels, n_latent
            h = self.decoder(gamma_reshape_, enum_label.to(z.device), batch_index_sc_input)
            px_est += proportion_modes[mode, ...].unsqueeze(-1) * self.px_decoder(h).reshape(
                (m, self.n_labels, -1)
            )  # (minibatch, n_labels, n_genes)

        # add the dummy cell type
        eps = eps.repeat((m, 1)).view(
            m, 1, -1
        )  # (M, 1, n_genes) <- this is the dummy cell type

        # account for gene specific bias and add noise
        r_hat = torch.cat(
            [beta.unsqueeze(0).unsqueeze(1) * px_est, eps], dim=1
        )  # M, n_labels + 1, n_genes
        # now combine them for convolution
        px_scale = torch.sum(v_ind.unsqueeze(2) * r_hat, dim=1)  # batch_size, n_genes
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

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        n_obs: int = 1.0,
        weighting_cross_entropy: float = 1e-6,
    ):
        """Compute the loss."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        gamma = generative_outputs["gamma"]
        v = generative_outputs["v"]

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)

        # eta prior likelihood
        mean = torch.zeros_like(self.eta)
        scale = torch.ones_like(self.eta)
        glo_neg_log_likelihood_prior = (
            -self.eta_reg * Normal(mean, scale).log_prob(self.eta).sum()
        )
        glo_neg_log_likelihood_prior += self.beta_reg * torch.var(self.beta)

        v_sparsity_loss = self.l1_reg * torch.abs(v).mean(1)

        # gamma prior likelihood
        if self.mean_vprior is None:
            # isotropic normal prior
            mean = torch.zeros_like(gamma)
            scale = torch.ones_like(gamma)
            neg_log_likelihood_prior = -Normal(mean, scale).log_prob(gamma).sum(2).sum(1)
        else:
            # vampprior
            # gamma is of shape minibatch_size, 1, n_latent, n_labels
            gamma = gamma.permute(3, 0, 2, 1) # minibatch_size, 1, n_labels, n_latent
            mean_vprior = torch.transpose(self.mean_vprior, 0, 1).unsqueeze(
                0
            )  # 1, p, n_labels, n_latent
            var_vprior = torch.transpose(self.var_vprior, 0, 1).unsqueeze(
                0
            )  # 1, p, n_labels, n_latent
            mp_vprior = torch.transpose(self.mp_vprior, 0, 1)  # p, n_labels
            if self.prior_mode == "mog":
                proportion_modes_logits = generative_outputs["proportion_modes_logits"]
                proportion_modes = generative_outputs["proportion_modes"]
                pre_lse = (
                    Normal(mean_vprior, torch.sqrt(var_vprior) + 1e-4).log_prob(gamma).sum(3)
                ) + torch.log(proportion_modes).permute(1, 0, 2)  # minibatch, p, n_labels
                log_likelihood_prior = torch.logsumexp(pre_lse, 1)  # minibatch, n_labels
                neg_log_likelihood_prior = -log_likelihood_prior.sum(1)  # minibatch

                neg_log_likelihood_prior += weighting_cross_entropy * torch.nn.functional.cross_entropy(
                    proportion_modes_logits.permute(1, 0, 2), mp_vprior.repeat(x.shape[0], 1, 1), reduction='none').sum(1)
            else:
                pre_lse = (
                    Normal(mean_vprior, torch.sqrt(var_vprior) + 1e-4).log_prob(gamma).sum(3)
                ) + torch.log(mp_vprior)  # minibatch, p, n_labels
                log_likelihood_prior = torch.logsumexp(pre_lse, 1)  # minibatch, n_labels
                neg_log_likelihood_prior = -log_likelihood_prior.sum(1)  # minibatch

        if self.n_latent_amortization is not None:
            neg_log_likelihood_prior += kl(
                inference_outputs["qz"],
                Normal(torch.zeros([self.n_latent_amortization], device=x.device), torch.ones([self.n_latent_amortization], device=x.device))
            ).sum(dim=-1)

        # High v_sparsity_loss is detrimental early in training, scaling by kl_weight to increase over training epochs.
        loss = n_obs * (
            torch.mean(
                reconst_loss + kl_weight * (neg_log_likelihood_prior + v_sparsity_loss)
            )
            + glo_neg_log_likelihood_prior
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

 
    @torch.inference_mode()
    @auto_move_data
    def get_ct_specific_expression(
        self, x: torch.Tensor = None, ind_x: torch.Tensor = None, y: int = None,
        batch_index: torch.Tensor = None, batch_index_sc: torch.Tensor = None
    ):
        """Returns cell type specific gene expression at the queried spots.

        Parameters
        ----------
        x
            tensor of data
        ind_x
            tensor of indices
        y
            integer for cell types
        batch_index_sc
            tensor of corresponding batch in single cell data for decoder
        """
        # cell-type specific gene expression, shape (minibatch, celltype, gene).
        beta = torch.exp(self.beta)  # n_genes
        y_torch = (y * torch.ones_like(ind_x)).ravel()
        # obtain the relevant gammas
        if self.amortization in ["both", "latent"]:
            x_ = torch.log(1 + x)
            z_amortization = self.amortization_network(x_, batch_index)[:, :self.n_latent_amortization]
            if self.prior_mode == "mog":
                gamma_ = self.gamma_encoder(z_amortization)
                proportion_modes_logits = torch.transpose(
                    gamma_[:, -self.p*self.n_labels:], 0, 1).reshape(
                    (self.p, 1, -1)
                ).transpose(1, 2)
                proportion_modes = torch.nn.functional.softmax(proportion_modes_logits, dim=0)
                # shape (p, n_labels, minibatch_size)
                gamma_ind = torch.transpose(
                    gamma_[:, :self.p*self.n_labels*self.n_latent], 0, 1).reshape(
                        (self.p, self.n_latent, self.n_labels, -1)
                )
            else:
                gamma_ind = torch.transpose(
                    self.gamma_encoder(z_amortization), 0, 1).reshape(
                        (1, self.n_latent, self.n_labels, -1)
                )
                proportion_modes = torch.ones((1, self.n_labels), device=x.device)
        else:
            gamma_ind = self.gamma[:, :, ind_x].unsqueeze(0)  # 1, n_latent, n_labels, minibatch_size
            proportion_modes = torch.ones((1, self.n_labels), device=x.device)
        gamma_ind = gamma_ind[:, :, y, :]
        proportion_modes = proportion_modes[:, y]

        px_est = torch.zeros((x.shape[0], self.n_genes), device=x.device)
        for mode in range(gamma_ind.shape[0]):
            # reshape and get gene expression value for all minibatch
            gamma_ind_ = torch.transpose(
                gamma_ind[mode, ...], 1, 0
            )  # minibatch_size, n_latent
            h = self.decoder(gamma_ind_, y_torch.unsqueeze(1), batch_index_sc.unsqueeze(1))
            px_est += proportion_modes[mode, ...].unsqueeze(-1) * self.px_decoder(h).reshape(
                (x.shape[0], -1)
            )  # (minibatch, n_genes)

        px_scale_ct = torch.exp(self.px_o).unsqueeze(0) * beta.unsqueeze(0) * px_est
        return px_scale_ct  # shape (minibatch, genes)

