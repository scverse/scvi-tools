import numpy as np
import torch

# import pdb
import torch.nn.functional as F

from scvi import _CONSTANTS
from scvi.compose import BaseModuleClass, LossRecorder, auto_move_data
from scvi.distributions import NegativeBinomial
from torch.distributions import Dirichlet, Normal


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
    **model_kwargs
        Additional kwargs
    """

    def __init__(
        self,
        n_genes: int,
        n_labels: int,
        rho: torch.Tensor,
        **model_kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = n_labels
        self.rho = rho

        self.register_buffer("cell_type_markers", rho)

        # perform all other initialization
        self.LOWER_BOUND = 1e-10
        self.THETA_LOWER_BOUND = 1e-20
        self.B = 10
        self.min_delta = 2
        self.dirichlet_concentration = [1e-2] * self.n_labels
        self.shrinkage = True

        # compute theta
        self.theta_logit = torch.nn.Parameter(torch.randn(self.n_labels))

        # compute delta (cell type specific overexpression parameter)
        self.delta_log = torch.nn.Parameter(
            torch.FloatTensor(self.n_genes, self.n_labels).uniform_(
                np.log(self.min_delta), 2
            )
        )

        # shrinkage prior on delta
        if self.shrinkage:
            self.delta_log_mean = torch.nn.Parameter(torch.Tensor(0))
            self.delta_log_variance = torch.nn.Parameter(torch.Tensor(1))

        self.log_a = torch.nn.Parameter(torch.zeros(self.B, dtype=torch.float64))

    def _get_inference_input(self, tensors):
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
        # x has shape (n, g)
        # n = 128
        # g = 100
        # c = 5
        delta = torch.exp(self.delta_log)  # (g, c)
        theta_log = F.log_softmax(self.theta_logit)  # (c)

        # compute mean of NegBin - shape (n_cells, n_genes, n_labels)
        s = x.sum(1, keepdim=True)  # (n, 1)
        base_mean = torch.log(s)
        base_mean_u = base_mean.unsqueeze(-1)  # (n, 1, 1)
        base_mean_e = base_mean_u.expand(
            s.shape[0], self.n_genes, self.n_labels
        )  # (n, g, c)

        delta_rho = delta * self.rho
        delta_rho_e = delta_rho.expand(
            s.shape[0], self.n_genes, self.n_labels
        )  # (n, g, c)
        log_mu_ngc = base_mean_e + delta_rho_e
        mu_ngc = torch.exp(log_mu_ngc)  # (n, g, c)

        # compute basis means for phi - shape (B)
        basis_means_fixed = np.linspace(torch.min(x), torch.max(x), self.B)
        basis_means = torch.tensor(basis_means_fixed)  # (B)

        # compute phi of NegBin - shape (n_cells, n_genes, n_labels)
        a = torch.exp(self.log_a)  # (B)
        a_e = a.expand(s.shape[0], self.n_genes, self.n_labels, self.B)
        b_init = 2 * ((basis_means_fixed[1] - basis_means_fixed[0]) ** 2)
        b = torch.exp(torch.ones(self.B) * (-np.log(b_init)))  # (B)
        b_e = b.expand(s.shape[0], self.n_genes, self.n_labels, self.B)
        mu_ngc_u = mu_ngc.unsqueeze(-1)
        mu_ngcb = mu_ngc_u.expand(
            s.shape[0], self.n_genes, self.n_labels, self.B
        )  # (n, g, c, B)
        basis_means_e = basis_means.expand(
            s.shape[0], self.n_genes, self.n_labels, self.B
        )  # (n, g, c, B)
        phi = (  # (n, g, c)
            torch.sum(a_e * torch.exp(-b_e * torch.square(mu_ngcb - basis_means_e)), 3)
            + self.LOWER_BOUND
        )

        # compute gamma
        nb_pdf = NegativeBinomial(probs=mu_ngc, total_count=phi)
        y_ = x.unsqueeze(-1).expand(s.shape[0], self.n_genes, self.n_labels)
        y_log_prob_raw = nb_pdf.log_prob(y_)  # (n, g, c)
        theta_log_e = theta_log.expand(s.shape[0], self.n_labels)
        p_y_unorm = torch.sum(y_log_prob_raw, 1) + theta_log_e  # (n, c)
        p_y_norm = torch.logsumexp(p_y_unorm, 1)
        p_y_norm_e = p_y_norm.unsqueeze(-1).expand(s.shape[0], self.n_labels)
        gamma = torch.exp(p_y_unorm - p_y_norm_e)  # (n, c)

        return dict(
            mu=mu_ngc,
            phi=phi,
            gamma=gamma,
            p_y_unorm=p_y_unorm,
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
        p_y_unorm = generative_outputs["p_y_unorm"]
        gamma = generative_outputs["gamma"]

        # compute Q
        # take mean of number of cells and multiply by n_obs (instead of summing n)
        q_per_cell = torch.sum(gamma * p_y_unorm, 1)

        # third term is log prob of prior terms in Q
        theta_log = F.log_softmax(self.theta_logit)
        theta_log_prior = Dirichlet(torch.tensor(self.dirichlet_concentration))
        theta_log_prob = -theta_log_prior.log_prob(
            torch.exp(theta_log) + self.THETA_LOWER_BOUND
        )
        delta_log_prior = Normal(self.delta_log_mean, self.delta_log_variance)
        summed_delta_log = torch.sum(self.delta_log)
        delta_log_prob = -torch.sum(delta_log_prior.log_prob(summed_delta_log))
        prior_log_prob = theta_log_prob
        if self.shrinkage:
            prior_log_prob += delta_log_prob

        loss = torch.mean(q_per_cell) * n_obs + prior_log_prob

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
