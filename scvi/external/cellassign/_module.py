from typing import Iterable, Optional

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Dirichlet, Normal

from scvi import _CONSTANTS
from scvi.compose import BaseModuleClass, LossRecorder, auto_move_data
from scvi.distributions import NegativeBinomial
from scvi.modules._utils import one_hot

LOWER_BOUND = 1e-10
THETA_LOWER_BOUND = 1e-20
B = 10


class CellAssignModule(BaseModuleClass):
    """
    Model for CellAssign.

    Parameters
    ----------
    n_genes
        Number of input genes
    n_labels
        Number of input cell types
    rho
        Binary matrix of cell type markers
    n_batch
        Number of batches, if 0, no batch correction is performed.
    n_cats_per_cov
        Number of categories for each extra categorical covariate
    n_continuous_cov
        Number of continuous covarites
    """

    def __init__(
        self,
        n_genes: int,
        rho: torch.Tensor,
        n_batch: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        n_continuous_cov: int = 0,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = rho.shape[1]
        self.n_batch = n_batch
        self.n_cats_per_cov = n_cats_per_cov
        self.n_continuous_cov = n_continuous_cov

        # this would be P in their code
        # if it's zero, let's ignore
        # though we do need the base expression intercept either way
        design_matrix_col_dim = n_batch + n_continuous_cov
        design_matrix_col_dim += 0 if n_cats_per_cov is None else sum(n_cats_per_cov)
        self.register_buffer("rho", rho)

        # perform all other initialization
        self.min_delta = 2
        dirichlet_concentration = torch.tensor([1e-2] * self.n_labels)
        self.register_buffer("dirichlet_concentration", dirichlet_concentration)
        self.shrinkage = True
        self.b_g_0 = torch.nn.Parameter(torch.randn(n_genes))

        # compute theta
        self.theta_logit = torch.nn.Parameter(torch.randn(self.n_labels))

        # compute delta (cell type specific overexpression parameter)
        self.delta_log_unclamped = torch.nn.Parameter(
            torch.FloatTensor(self.n_genes, self.n_labels).uniform_(-2, 2)
        )

        # shrinkage prior on delta
        if self.shrinkage:
            self.delta_log_mean = torch.nn.Parameter(
                torch.zeros(
                    1,
                )
            )
            self.delta_log_variance = torch.nn.Parameter(
                torch.ones(
                    1,
                )
            )

        self.log_a = torch.nn.Parameter(torch.zeros(B, dtype=torch.float64))

    def _get_inference_input(self, tensors):
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[_CONSTANTS.X_KEY]
        size_factor = tensors["_size_factor"]

        to_cat = []
        if self.n_batch > 0:
            to_cat.append(one_hot(tensors[_CONSTANTS.BATCH_KEY], self.n_batch))

        cont_key = _CONSTANTS.CONT_COVS_KEY
        if cont_key in tensors.keys():
            to_cat.append(tensors[cont_key])

        cat_key = _CONSTANTS.CAT_COVS_KEY
        if cat_key in tensors.keys():
            for cat_input, n_cat in zip(
                torch.split(tensors[cat_key], 1, dim=1), self.n_cats_per_cov
            ):
                to_cat.append(one_hot(cat_input, n_cat))

        design_matrix = torch.cat(to_cat) if len(to_cat) > 0 else None

        input_dict = dict(x=x, size_factor=size_factor, design_matrix=design_matrix)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, size_factor, design_matrix=None):
        self.delta_log = torch.clamp(
            self.delta_log_unclamped, min=np.log(self.min_delta)
        )
        # x has shape (n, g)
        delta = torch.exp(self.delta_log)  # (g, c)
        theta_log = F.log_softmax(self.theta_logit)  # (c)

        # compute mean of NegBin - shape (n_cells, n_genes, n_labels)
        s = size_factor
        n_cells = s.shape[0]
        base_mean = torch.log(s)  # (n, 1)
        base_mean_u = base_mean.unsqueeze(-1)  # (n, 1, 1)
        base_mean_e = base_mean_u.expand(
            n_cells, self.n_genes, self.n_labels
        )  # (n, g, c)

        # compute beta (covariate coefficent)
        # design_matrix has shape (n,p)
        if design_matrix is not None:
            col_means = torch.mean(x, 0)  # (g)
            col_means_mu, col_means_std = torch.std_mean(col_means)
            beta_0_init = torch.div(torch.sub(col_means, col_means_mu), col_means_std)
            beta_init = torch.zeros(
                [self.n_genes, design_matrix.shape[1] - 1]
            )  # (g, p-1)
            self.beta = torch.nn.Parameter(
                torch.cat((beta_0_init.unsqueeze(-1), beta_init), 1)
            )  # (g, p)
            covariates = torch.einsum("np,gp->gn", design_matrix, self.beta)  # (g, n)
            covariates_u = torch.transpose(covariates, 0, 1).unsqueeze(-1)  # (n, g, 1)
            covariates_e = covariates_u.expand(n_cells, self.n_genes, self.n_labels)
            base_mean_e = base_mean_e + covariates_e

        # base gene expression
        b_g_0_e = self.b_g_0.unsqueeze(-1).expand(n_cells, self.n_genes, self.n_labels)

        delta_rho = delta * self.rho
        delta_rho_e = delta_rho.expand(
            n_cells, self.n_genes, self.n_labels
        )  # (n, g, c)
        log_mu_ngc = base_mean_e + delta_rho_e + b_g_0_e
        mu_ngc = torch.exp(log_mu_ngc)  # (n, g, c)

        # compute basis means for phi - shape (B)
        basis_means = torch.linspace(
            torch.min(x), torch.max(x), B, device=x.device
        )  # (B)

        # compute phi of NegBin - shape (n_cells, n_genes, n_labels)
        a = torch.exp(self.log_a)  # (B)
        a_e = a.expand(n_cells, self.n_genes, self.n_labels, B)
        b_init = 2 * ((basis_means[1] - basis_means[0]) ** 2)
        b = torch.exp(torch.ones(B, device=x.device) * (-torch.log(b_init)))  # (B)
        b_e = b.expand(n_cells, self.n_genes, self.n_labels, B)
        mu_ngc_u = mu_ngc.unsqueeze(-1)
        mu_ngcb = mu_ngc_u.expand(
            n_cells, self.n_genes, self.n_labels, B
        )  # (n, g, c, B)
        basis_means_e = basis_means.expand(
            n_cells, self.n_genes, self.n_labels, B
        )  # (n, g, c, B)
        phi = (  # (n, g, c)
            torch.sum(a_e * torch.exp(-b_e * torch.square(mu_ngcb - basis_means_e)), 3)
            + LOWER_BOUND
        )

        # compute gamma
        nb_pdf = NegativeBinomial(mu=mu_ngc, theta=phi)
        y_ = x.unsqueeze(-1).expand(n_cells, self.n_genes, self.n_labels)
        y_log_prob_raw = nb_pdf.log_prob(y_)  # (n, g, c)
        theta_log_e = theta_log.expand(n_cells, self.n_labels)
        p_y_c = torch.sum(y_log_prob_raw, 1) + theta_log_e  # (n, c)
        normalizer_over_c = torch.logsumexp(p_y_c, 1)
        normalizer_over_c_e = normalizer_over_c.unsqueeze(-1).expand(
            n_cells, self.n_labels
        )
        gamma = torch.exp(p_y_c - normalizer_over_c_e)  # (n, c)

        return dict(
            mu=mu_ngc,
            phi=phi,
            gamma=gamma,
            p_y_c=p_y_c,
            s=s,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        n_obs: int = 1.0,
    ):
        # generative_outputs is a dict of the return value from `generative(...)`
        # assume that `n_obs` is the number of training data points
        p_y_c = generative_outputs["p_y_c"]
        gamma = generative_outputs["gamma"]

        # compute Q
        # take mean of number of cells and multiply by n_obs (instead of summing n)
        q_per_cell = torch.sum(gamma * -p_y_c, 1)

        # third term is log prob of prior terms in Q
        theta_log = F.log_softmax(self.theta_logit)
        theta_log_prior = Dirichlet(self.dirichlet_concentration)
        theta_log_prob = -theta_log_prior.log_prob(
            torch.exp(theta_log) + THETA_LOWER_BOUND
        )
        prior_log_prob = theta_log_prob
        if self.shrinkage:
            delta_log_prior = Normal(self.delta_log_mean, self.delta_log_variance)
            delta_log_prob = torch.masked_select(
                delta_log_prior.log_prob(self.delta_log), (self.rho > 0)
            )
            prior_log_prob += -torch.sum(delta_log_prob)

        loss = torch.mean(q_per_cell) * n_obs + prior_log_prob / n_obs

        return LossRecorder(
            loss, q_per_cell, torch.zeros_like(q_per_cell), prior_log_prob
        )

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for CellAssign")
