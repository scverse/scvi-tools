from typing import Tuple

import numpy as np
import torch
from torch.distributions import NegativeBinomial, Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.compose import BaseModuleClass, LossRecorder, auto_move_data


class RNADeconv(BaseModuleClass):
    """
    Model of single-cell RNA-sequencing data for deconvolution of spatial transriptomics.

    Reimplementation of the ScModel module of Stereoscope [Andersson20]_:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_genes
        Number of input genes
    n_labels
        Number of input cell types
    **model_kwargs
        Additional kwargs
    """

    def __init__(
        self,
        n_genes: int,
        n_labels: int,
        **model_kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = n_labels

        # logit param for negative binomial
        self.px_o = torch.nn.Parameter(torch.randn(self.n_genes))
        self.W = torch.nn.Parameter(
            torch.randn(self.n_genes, self.n_labels)
        )  # n_genes, n_cell types

        if "ct_weight" in model_kwargs:
            ct_weight = torch.tensor(model_kwargs["ct_prop"], dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)

    @torch.no_grad()
    def get_params(self) -> Tuple[np.ndarray]:
        """
        Returns the parameters for feeding into the spatial data
        -------
        type
            list of tensor
        """
        return self.W.cpu().numpy(), self.px_o.cpu().numpy()

    def _get_inference_input(self, tensors):
        # we perform MAP here, so there is nothing to infer
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[_CONSTANTS.X_KEY]
        y = tensors[_CONSTANTS.LABELS_KEY]

        input_dict = dict(x=x, y=y)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, y):
        """Simply build the negative binomial parameters for every cell in the minibatch."""
        px_scale = torch.nn.functional.softplus(self.W)[
            :, y.long()[:, 0]
        ].T  # cells per gene
        library = torch.sum(x, dim=1, keepdim=True)
        px_rate = library * px_scale
        scaling_factor = self.ct_weight[y.long()[:, 0]]

        return dict(
            px_scale=px_scale,
            px_o=self.px_o,
            px_rate=px_rate,
            library=library,
            scaling_factor=scaling_factor,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        x = tensors[_CONSTANTS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        scaling_factor = generative_outputs["scaling_factor"]

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
        loss = torch.mean(scaling_factor * reconst_loss)

        return LossRecorder(loss, reconst_loss, torch.zeros((1,)), 0.0)

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for Stereoscope")


class SpatialDeconv(BaseModuleClass):
    """
    Model of single-cell RNA-sequencing data for deconvolution of spatial transriptomics.

    Reimplementation of the STModel module of Stereoscope [Andersson20]_:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_spots
        Number of input spots
    sc_params
        Tuple of ndarray of shapes [(n_genes, n_labels), (n_genes)] containing the dictionnary and log dispersion parameters
    prior_weight
        Whether to sample the minibatch by the number of total observations or the monibatch size
    """

    def __init__(
        self,
        n_spots: int,
        sc_params: Tuple[np.ndarray],
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
    ):
        super().__init__()
        # unpack and copy parameters
        w, px_o = sc_params
        self.register_buffer("W", torch.tensor(w))
        self.register_buffer("px_o", torch.tensor(px_o))

        # setup constants
        self.n_spots = n_spots
        self.n_genes, self.n_labels = self.W.shape
        self.prior_weight = prior_weight

        # noise from data
        self.eta = torch.nn.Parameter(torch.randn(self.n_genes))
        # factor loadings
        self.V = torch.nn.Parameter(torch.randn(self.n_labels + 1, self.n_spots))
        # additive gene bias
        self.beta = torch.nn.Parameter(0.01 * torch.randn(self.n_genes))

    @torch.no_grad()
    def get_proportions(self, keep_noise=False) -> np.ndarray:
        """Returns the loadings."""
        # get estimated unadjusted proportions
        res = (
            torch.nn.functional.softplus(self.V).cpu().numpy().T
        )  # n_spots, n_labels + 1
        # remove dummy cell type proportion values
        if not keep_noise:
            res = res[:, :-1]
        # normalize to obtain adjusted proportions
        res = res / res.sum(axis=1).reshape(-1, 1)
        return res

    def _get_inference_input(self, tensors):
        # we perform MAP here, so there is nothing to infer
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[_CONSTANTS.X_KEY]
        ind_x = tensors["ind_x"]

        input_dict = dict(x=x, ind_x=ind_x)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, ind_x):
        """Build the deconvolution model for every cell in the minibatch."""

        beta = torch.nn.functional.softplus(self.beta)  # n_genes
        v = torch.nn.functional.softplus(self.V)  # n_labels + 1, n_spots
        w = torch.nn.functional.softplus(self.W)  # n_genes, n_labels
        eps = torch.nn.functional.softplus(self.eta)  # n_genes

        # account for gene specific bias and add noise
        r_hat = torch.cat(
            [beta.unsqueeze(1) * w, eps.unsqueeze(1)], dim=1
        )  # n_genes, n_labels + 1
        # subsample observations
        v_ind = v[:, ind_x.long()[:, 0]]  # labels + 1, batch_size
        px_rate = torch.transpose(
            torch.matmul(r_hat, v_ind), 0, 1
        )  # batch_size, n_genes

        return dict(px_o=self.px_o, px_rate=px_rate, eta=self.eta)

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

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
        # prior likelihood
        mean = torch.zeros_like(self.eta)
        scale = torch.ones_like(self.eta)
        neg_log_likelihood_prior = -Normal(mean, scale).log_prob(self.eta).sum()

        if self.prior_weight == "n_obs":
            # the correct way to reweight observations while performing stochastic optimization
            loss = n_obs * torch.mean(reconst_loss) + neg_log_likelihood_prior
        else:
            # the original way it is done in Stereoscope; we use this option to show reproducibility of their codebase
            loss = torch.sum(reconst_loss) + neg_log_likelihood_prior
        return LossRecorder(
            loss, reconst_loss, torch.zeros((1,)), neg_log_likelihood_prior
        )

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for Stereoscope")
