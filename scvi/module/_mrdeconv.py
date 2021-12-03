from collections import OrderedDict

import numpy as np
import torch
from torch.distributions import Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import FCLayers


def identity(x):
    return x


class MRDeconv(BaseModuleClass):
    """
    Model for multi-resolution deconvolution of spatial transriptomics.

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
    amortization
        which of the latent variables to amortize inference over (gamma, proportions, both or none)
    """

    def __init__(
        self,
        n_spots: int,
        n_labels: int,
        n_hidden: int,
        n_layers: int,
        n_latent: int,
        n_genes: int,
        decoder_state_dict: OrderedDict,
        px_decoder_state_dict: OrderedDict,
        px_r: np.ndarray,
        mean_vprior: np.ndarray = None,
        var_vprior: np.ndarray = None,
        amortization: Literal["none", "latent", "proportion", "both"] = "latent",
    ):
        super().__init__()
        self.n_spots = n_spots
        self.n_labels = n_labels
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.n_genes = n_genes
        self.amortization = amortization
        # unpack and copy parameters
        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=[n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            use_layer_norm=True,
            use_batch_norm=False,
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
        else:
            self.mean_vprior = None
            self.var_vprior = None
        # noise from data
        self.eta = torch.nn.Parameter(torch.randn(self.n_genes))
        # additive gene bias
        self.beta = torch.nn.Parameter(0.01 * torch.randn(self.n_genes))

        # create additional neural nets for amortization
        # within cell_type factor loadings
        self.gamma_encoder = torch.nn.Sequential(
            FCLayers(
                n_in=self.n_genes,
                n_out=n_hidden,
                n_cat_list=None,
                n_layers=2,
                n_hidden=n_hidden,
                dropout_rate=0.1,
            ),
            torch.nn.Linear(n_hidden, n_latent * n_labels),
        )
        # cell type loadings
        self.V_encoder = FCLayers(
            n_in=self.n_genes,
            n_out=self.n_labels + 1,
            n_layers=2,
            n_hidden=n_hidden,
            dropout_rate=0.1,
        )

    def _get_inference_input(self, tensors):
        # we perform MAP here, so we just need to subsample the variables
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[_CONSTANTS.X_KEY]
        ind_x = tensors["ind_x"].long()

        input_dict = dict(x=x, ind_x=ind_x)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, ind_x):
        """Build the deconvolution model for every cell in the minibatch."""
        m = x.shape[0]
        library = torch.sum(x, dim=1, keepdim=True)
        # setup all non-linearities
        beta = torch.nn.functional.softplus(self.beta)  # n_genes
        eps = torch.nn.functional.softplus(self.eta)  # n_genes
        x_ = torch.log(1 + x)
        # subsample parameters

        if self.amortization in ["both", "latent"]:
            gamma_ind = torch.transpose(self.gamma_encoder(x_), 0, 1).reshape(
                (self.n_latent, self.n_labels, -1)
            )
        else:
            gamma_ind = self.gamma[
                :, :, ind_x[:, 0]
            ]  # n_latent, n_labels, minibatch_size

        if self.amortization in ["both", "proportion"]:
            v_ind = self.V_encoder(x_)
        else:
            v_ind = self.V[:, ind_x[:, 0]].T  # minibatch_size, labels + 1
        v_ind = torch.nn.functional.softplus(v_ind)

        # reshape and get gene expression value for all minibatch
        gamma_ind = torch.transpose(
            gamma_ind, 2, 0
        )  # minibatch_size, n_labels, n_latent
        gamma_reshape = gamma_ind.reshape(
            (-1, self.n_latent)
        )  # minibatch_size * n_labels, n_latent
        enum_label = (
            torch.arange(0, self.n_labels).repeat((m)).view((-1, 1))
        )  # minibatch_size * n_labels, 1
        h = self.decoder(gamma_reshape, enum_label.to(x.device))
        px_rate = self.px_decoder(h).reshape(
            (m, self.n_labels, -1)
        )  # (minibatch, n_labels, n_genes)

        # add the dummy cell type
        eps = eps.repeat((m, 1)).view(
            m, 1, -1
        )  # (M, 1, n_genes) <- this is the dummy cell type

        # account for gene specific bias and add noise
        r_hat = torch.cat(
            [beta.unsqueeze(0).unsqueeze(1) * px_rate, eps], dim=1
        )  # M, n_labels + 1, n_genes
        # now combine them for convolution
        px_scale = torch.sum(v_ind.unsqueeze(2) * r_hat, dim=1)  # batch_size, n_genes
        px_rate = library * px_scale

        return dict(
            px_o=self.px_o, px_rate=px_rate, px_scale=px_scale, gamma=gamma_ind, v=v_ind
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        n_obs: int = 1.0,
    ):
        x = tensors[_CONSTANTS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        gamma = generative_outputs["gamma"]

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)

        # eta prior likelihood
        mean = torch.zeros_like(self.eta)
        scale = torch.ones_like(self.eta)
        glo_neg_log_likelihood_prior = -Normal(mean, scale).log_prob(self.eta).sum()
        glo_neg_log_likelihood_prior += torch.var(self.beta)

        # gamma prior likelihood
        if self.mean_vprior is None:
            # isotropic normal prior
            mean = torch.zeros_like(gamma)
            scale = torch.ones_like(gamma)
            neg_log_likelihood_prior = (
                -Normal(mean, scale).log_prob(gamma).sum(2).sum(1)
            )
        else:
            # vampprior
            # gamma is of shape n_latent, n_labels, minibatch_size
            gamma = gamma.unsqueeze(1)  # minibatch_size, 1, n_labels, n_latent
            mean_vprior = torch.transpose(self.mean_vprior, 0, 1).unsqueeze(
                0
            )  # 1, p, n_labels, n_latent
            var_vprior = torch.transpose(self.var_vprior, 0, 1).unsqueeze(
                0
            )  # 1, p, n_labels, n_latent
            pre_lse = (
                Normal(mean_vprior, torch.sqrt(var_vprior)).log_prob(gamma).sum(-1)
            )  # minibatch, p, n_labels
            log_likelihood_prior = torch.logsumexp(pre_lse, 1) - np.log(
                self.p
            )  # minibatch, n_labels
            neg_log_likelihood_prior = -log_likelihood_prior.sum(1)  # minibatch
            # mean_vprior is of shape n_labels, p, n_latent

        loss = (
            n_obs * torch.mean(reconst_loss + kl_weight * neg_log_likelihood_prior)
            + glo_neg_log_likelihood_prior
        )

        return LossRecorder(
            loss, reconst_loss, neg_log_likelihood_prior, glo_neg_log_likelihood_prior
        )

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for DestVI")

    @torch.no_grad()
    @auto_move_data
    def get_proportions(self, x=None, keep_noise=False) -> np.ndarray:
        """Returns the loadings."""
        if self.amortization in ["both", "proportion"]:
            # get estimated unadjusted proportions
            x_ = torch.log(1 + x)
            res = torch.nn.functional.softplus(self.V_encoder(x_))
        else:
            res = (
                torch.nn.functional.softplus(self.V).cpu().numpy().T
            )  # n_spots, n_labels + 1
        # remove dummy cell type proportion values
        if not keep_noise:
            res = res[:, :-1]
        # normalize to obtain adjusted proportions
        res = res / res.sum(axis=1).reshape(-1, 1)
        return res

    @torch.no_grad()
    @auto_move_data
    def get_gamma(self, x: torch.Tensor = None) -> torch.Tensor:
        """
        Returns the loadings.

        Returns
        -------
        type
            tensor
        """
        # get estimated unadjusted proportions
        if self.amortization in ["latent", "both"]:
            x_ = torch.log(1 + x)
            gamma = self.gamma_encoder(x_)
            return torch.transpose(gamma, 0, 1).reshape(
                (self.n_latent, self.n_labels, -1)
            )  # n_latent, n_labels, minibatch
        else:
            return self.gamma.cpu().numpy()  # (n_latent, n_labels, n_spots)

    @torch.no_grad()
    @auto_move_data
    def get_ct_specific_expression(
        self, x: torch.Tensor = None, ind_x: torch.Tensor = None, y: int = None
    ):
        """
        Returns cell type specific gene expression at the queried spots.

        Parameters
        ----------
        x
            tensor of data
        ind_x
            tensor of indices
        y
            integer for cell types
        """
        # cell-type specific gene expression, shape (minibatch, celltype, gene).
        beta = torch.nn.functional.softplus(self.beta)  # n_genes
        y_torch = y * torch.ones_like(ind_x)
        # obtain the relevant gammas
        if self.amortization in ["both", "latent"]:
            x_ = torch.log(1 + x)
            gamma_ind = torch.transpose(self.gamma_encoder(x_), 0, 1).reshape(
                (self.n_latent, self.n_labels, -1)
            )
        else:
            gamma_ind = self.gamma[
                :, :, ind_x[:, 0]
            ]  # n_latent, n_labels, minibatch_size

        # calculate cell type specific expression
        gamma_select = gamma_ind[
            :, y_torch[:, 0], torch.arange(ind_x.shape[0])
        ].T  # minibatch_size, n_latent
        h = self.decoder(gamma_select, y_torch)
        px_scale = self.px_decoder(h)  # (minibatch, n_genes)
        px_ct = torch.exp(self.px_o).unsqueeze(0) * beta.unsqueeze(0) * px_scale
        return px_ct  # shape (minibatch, genes)
