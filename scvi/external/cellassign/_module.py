import numpy as np
import torch
import pdb

from scvi import _CONSTANTS
from scvi.compose import AbstractVAE, SCVILoss, auto_move_data
from scvi.distributions import NegativeBinomial
from torch.distributions import Dirichlet, Normal


class CellAssignModule(AbstractVAE):
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

        # self.register_buffer("rho", rho)

        # perform all other initialization
        self.LOWER_BOUND = 1e-10
        self.THETA_LOWER_BOUND = 1e-20
        self.B = 10
        self.min_delta = 2
        self.dirichlet_concentration = [1e-2] * self.n_labels
        self.shrinkage = True

        # compute theta
        self.theta_logit = torch.nn.Parameter(
            torch.randn(self.n_labels, dtype=torch.float64)
        )

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
        # n = 128
        # g = 100
        # c = 5
        pdb.set_trace()
        delta = torch.exp(self.delta_log)  # (100, 5)
        log_softmax = torch.nn.LogSoftmax()
        theta_log = log_softmax(self.theta_logit)  # (5)

        # compute mean of NegBin - shape (n_cells, n_genes, n_labels)
        s = x.sum(1, keepdim=True)  # (128, 1)
        base_mean = torch.log(s)
        base_mean_u = base_mean.unsqueeze(-1)  # (128, 1, 1)
        base_mean_e = base_mean_u.expand(
            s.shape[0], self.n_genes, self.n_labels
        )  # (128, 100, 5)

        delta_rho = delta * self.rho
        delta_rho_e = delta_rho.expand(
            s.shape[0], self.n_genes, self.n_labels
        )  # (128, 100, 5)

        log_mu_ngc = base_mean_e + delta_rho_e
        mu_ngc = torch.exp(log_mu_ngc)  # (128, 100, 5)

        # compute basis means for phi - shape (B)
        basis_means_fixed = np.linspace(torch.min(x), torch.max(x), self.B)
        basis_means = torch.tensor(basis_means_fixed)  # (10)

        # compute phi of NegBin - shape (n_cells, n_genes, n_labels)
        a = torch.exp(self.log_a)  # (10)
        a_e = a.expand(s.shape[0], self.n_genes, self.n_labels, self.B)
        b_init = 2 * ((basis_means_fixed[1] - basis_means_fixed[0]) ** 2)
        b = torch.exp(torch.ones(self.B) * (-np.log(b_init)))  # (10)
        b_e = b.expand(s.shape[0], self.n_genes, self.n_labels, self.B)
        mu_ngc_u = mu_ngc.unsqueeze(-1)
        mu_ngcb = mu_ngc_u.expand(
            s.shape[0], self.n_genes, self.n_labels, self.B
        )  # (128, 100, 5, 10)
        basis_means_e = basis_means.expand(
            s.shape[0], self.n_genes, self.n_labels, self.B
        )  # (128, 100, 5, 10)
        phi = (  # (128, 100, 5)
            torch.sum(a_e * torch.exp(-b_e * torch.square(mu_ngcb - basis_means_e)), 3)
            + self.LOWER_BOUND
        )

        # compute gamma
        p = mu_ngc / (mu_ngc + phi)
        nb_pdf = NegativeBinomial(probs=p, total_count=phi)
        y_tensor_list = []
        for _ in range(self.n_labels):
            y_tensor_list += [x]
        y_ = torch.stack(y_tensor_list, axis=2)
        y_log_prob_raw = nb_pdf.log_prob(y_)  # (128, 100, 5)
        theta_log_e = theta_log.expand(128, 5)
        p_y_unorm = torch.sum(y_log_prob_raw, 1) + theta_log_e  # (128, 5)
        p_y_norm = torch.logsumexp(p_y_unorm, 1)
        p_y_norm_e = p_y_norm.unsqueeze(-1).expand(128, 5)
        gamma = torch.exp(p_y_unorm - p_y_norm_e)  # (128, 5)

        return dict(
            mu=mu_ngc,
            phi=phi,
            gamma=gamma,
            p_y_norm=p_y_norm,
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
        p_y_norm = generative_outputs["p_y_norm"]

        # compute Q
        # gamma_fixed = torch.empty((None, self.n_labels))
        # Q = -torch.einsum("nc,cn->", gamma_fixed, p_y_unorm)
        Q = p_y_norm

        # second term in SCVILoss is Q per cell without prior terms. shape is (n_cells,)
        loss = torch.sum(p_y_norm)

        # third term is log prob of prior terms in Q
        log_softmax = torch.nn.LogSoftmax()
        theta_log = log_softmax(self.theta_logit)
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

        return SCVILoss(loss, Q, prior_log_prob, 0.0)

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for CellAssign")
