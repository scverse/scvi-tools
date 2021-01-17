import numpy as np
import torch

from scvi import _CONSTANTS
from scvi.compose import AbstractVAE, SCVILoss, auto_move_data
from scvi.distributions import NegativeBinomial


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
        # n_cells: int,
        # n_covariates: int,
        **model_kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = n_labels
        self.rho = rho
        # self.n_cells = n_cells
        # self.n_covariates = n_covariates

        self.register_buffer("rho", rho)

        # perform all other initialization
        self.LOWER_BOUND = 1e-10
        self.THETA_LOWER_BOUND = 1e-20
        self.B = 10
        self.min_delta = 2
        self.dirichlet_concentration = [1e-2] * self.n_labels
        # self.shrinkage = True
        # self.verbose = False
        # self.n_batches = 1
        # self.rel_tol_adam = 1e-4
        # self.rel_tol_em = 1e-4
        # self.max_iter_adam = 1e5
        # self.max_iter_em = 20
        # self.learning_rate = 1e-4
        # self.random_seed = None
        # self.threads = 0

        # compute theta
        theta_logit = torch.randn(self.n_labels, dtype=torch.float64)
        self.theta_log = torch.Parameter(torch.nn.LogSoftmax(theta_logit))

        # compute delta (cell type specific overexpression parameter)
        self.delta_log = torch.Parameter(
            torch.FloatTensor(self.n_genes, self.n_labels).uniform_(
                np.log(self.min_delta), 2
            )
        )

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

        # compute beta (covariate coefficent)
        # col_means = torch.mean(x, 0)
        # col_means_mu, col_means_std = torch.std_mean(col_means)
        # beta_0_init = torch.div(torch.sub(col_means, col_means_mu), col_means_std)
        # beta_init = torch.zeros([self.n_genes, self.n_covariates - 1], dtype=torch.float64)
        # beta = torch.cat((beta_0_init, beta_init), 1)

        delta = torch.exp(self.delta_log)

        # compute mean of NegBin
        s = x.sum(1, keepDim=True)
        base_mean = torch.transpose(torch.log(s))
        base_mean_list = []
        for _ in range(self.n_labels):
            base_mean_list += [base_mean]
        mu_ngc = torch.add(
            torch.stack(base_mean_list, 2), torch.multiply(delta, self.rho)
        )
        mu_cng = torch.transpose(mu_ngc, (2, 0, 1))
        mu_cngb = torch.tile(torch.expand_dims(mu_cng, axis=3), (1, 1, 1, self.B))
        mu_ngc = torch.transpose(mu_cng, (1, 2, 0))
        mu_ngc = torch.exp(mu_ngc)

        # compute basis means for phi
        basis_means_fixed = np.linspace(min(x), max(x), self.B)
        basis_means = torch.tensor(basis_means_fixed, dtype=torch.float64)

        # compute phi of NegBin
        a = torch.exp(torch.zeros(self.B, dtype=torch.float64))
        b_init = 2 * ((basis_means_fixed[1] - basis_means_fixed[0]) ** 2)
        b = torch.exp(torch.ones(self.B, dtype=torch.float64) * (-np.log(b_init)))
        phi_cng = (
            torch.reduce_sum(a * torch.exp(-b * torch.square(mu_cngb - basis_means)), 3)
            + self.LOWER_BOUND
        )
        phi = torch.transpose(phi_cng, (1, 2, 0))

        # compute gamma
        p = mu_ngc / (mu_ngc + phi)
        nb_pdf = NegativeBinomial(probs=p, total_count=phi)
        y_tensor_list = []
        for _ in range(self.n_labels):
            y_tensor_list += [x]
        y_ = torch.stack(y_tensor_list, axis=2)
        y_log_prob_raw = nb_pdf.log_prob(y_)
        y_log_prob = torch.transpose(y_log_prob_raw, (0, 2, 1))
        y_log_prob_sum = torch.reduce_sum(y_log_prob, 2) + self.theta_log
        p_y_on_c_unorm = torch.transpose(y_log_prob_sum, (1, 0))
        p_y_on_c_norm = torch.reshape(
            torch.reduce_logsumexp(p_y_on_c_unorm, 0), (1, -1)
        )
        gamma = torch.transpose(torch.exp(p_y_on_c_unorm - p_y_on_c_norm))

        return dict(
            mu=mu_ngc,
            phi=phi,
            gamma=gamma,
            p_y_on_c_unorm=p_y_on_c_unorm,
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
        # mu = generative_outputs['mu']
        # phi = generative_outputs['phi']
        # gamma = generative_outputs['gamma']
        p_y_on_c_unorm = generative_outputs["p_y_on_c_unorm"]

        # compute Q and put it in loss
        gamma_fixed = torch.empty((None, self.n_labels))
        Q = -torch.einsum("nc,cn->", gamma_fixed, p_y_on_c_unorm)

        # second term in SCVILoss is Q per cell without prior terms. shape is (n_cells,)
        loss = torch.reduce_sum(torch.reduce_logsumexp(p_y_on_c_unorm, 0))

        # third term is log prob of prior terms in Q
        # see the VAE class for how to compute NegBin log probability
        theta_log_prior = torch.Dirichlet(self.dirichlet_concentration)
        theta_log_prob = -theta_log_prior.log_prob(
            torch.exp(self.theta_log) + self.THETA_LOWER_BOUND
        )
        delta_log_prior = torch.Normal(0 * self.rho, 1)
        delta_log_prob = -torch.reduce_sum(delta_log_prior.log_prob(self.delta_log))
        prior_log_prob = theta_log_prob + delta_log_prob

        return SCVILoss(loss, Q, prior_log_prob, 0.0)

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for CellAssign")
