import numpy as np
import pandas as pd
import pyro
import pyro.distributions as dist
import torch

# import scvi
from pyro.distributions import constraints
from pyro.distributions.transforms import SoftplusTransform
from pyro.nn import PyroModule
from scipy.sparse import issparse
from torch.distributions import biject_to, transform_to

from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry
from scvi.distributions._negative_binomial import _convert_mean_disp_to_counts_logits
from scvi.nn import one_hot


@biject_to.register(constraints.positive)
@transform_to.register(constraints.positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()


# class NegativeBinomial(TorchDistributionMixin, ScVINegativeBinomial):
#    pass


class LocationModelLinearDependentWMultiExperimentModel(PyroModule):
    """
    Cell2location models the elements of :math:`D` as Negative Binomial distributed,
    given an unobserved gene expression level (rate) :math:`mu` and a gene- and batch-specific
    over-dispersion parameter :math:`\alpha_{e,g}` which accounts for unexplained variance:

    .. math::
        D_{s,g} \sim \mathtt{NB}(\mu_{s,g}, \alpha_{e,g})

    The expression level of genes :math:`\mu_{s,g}` in the mRNA count space is modelled
    as a linear function of expression signatures of reference cell types :math:`g_{f,g}`:

    .. math::
        \mu_{s,g} = m_{g} \left (\sum_{f} {w_{s,f} \: g_{f,g}} \right) + l_s + s_{e,g}

    Here, :math:`w_{s,f}` denotes regression weight of each reference signature :math:`f` at location :math:`s`,
      which can be interpreted as the expected number of cells at location :math:`s`
      that express reference signature :math:`f`;
    :math:`g_{f,g}` denotes the reference signatures of cell types :math:`f` of each gene :math:`g`,
      `cell_state_df` input ;
    :math:`m_{g}` denotes a gene-specific scaling parameter which adjusts for global differences in sensitivity
      between technologies (platform effect);
    :math:`l_{s}` and :math:`s_{e,g}` are additive components that account for gene- and location-specific shift,
      such as due to contaminating or free-floating RNA.

    To account for the similarity of location patterns across cell types, :math:`w_{s,f}` is modelled using
    another layer  of decomposition (factorization) using :math:`r={1, .., R}` groups of cell types,
    that can be interpreted as cellular compartments or tissue zones. Unless stated otherwise, R is set to 50.

    Corresponding graphical model can be found in supplementary methods:
    https://www.biorxiv.org/content/10.1101/2020.11.15.378125v1.supplementary-material

    Approximate Variational Inference is used to estimate the posterior distribution of all model parameters.

    Estimation of absolute cell abundance `w_{s,f}` is guided using informed prior on the number of cells
    (argument called `N_cells_per_location`). It is a tissue-level global estimate, which can be derived from histology
    images (H&E or DAPI), ideally paired to the spatial expression data or at least representing the same tissue type.
    This parameter can be estimated by manually counting nuclei in a 10-20 locations in the histology image
    (e.g. using 10X Loupe browser), and computing the average cell abundance.
    An appropriate setting of this prior is essential to inform the estimation of absolute cell type abundance values,
    however, the model is robust to a range of similar values.
    In settings where suitable histology images are not available, the size of capture regions relative to
    the expected size of cells can be used to estimate `N_cells_per_location`.
    """

    def __init__(
        self,
        n_obs,
        n_vars,
        n_factors,
        n_batch,
        cell_state_mat,
        batch_size=None,
        n_groups: int = 50,
        m_g_gene_level_prior={"mean": 1 / 2, "sd": 1 / 4},
        m_g_gene_level_var_prior={"mean_var_ratio": 1.0},
        N_cells_per_location=8.0,
        A_factors_per_location=7.0,
        Y_groups_per_location=7.0,
        N_cells_mean_var_ratio=1.0,
        alpha_g_phi_hyp_prior={"alpha": 9.0, "beta": 3.0},
        gene_add_alpha_hyp_prior={"alpha": 9.0, "beta": 3.0},
        gene_add_mean_hyp_prior={
            "alpha": 1.0,
            "beta": 100.0,
        },  # TODO initialise as average of empty locations
        w_sf_mean_var_ratio=5.0,
    ):

        super().__init__()

        self.n_obs = n_obs
        self.n_vars = n_vars
        self.n_factors = n_factors
        self.n_batch = n_batch
        self.batch_size = batch_size
        self.n_groups = n_groups

        for k in m_g_gene_level_var_prior.keys():
            m_g_gene_level_prior[k] = m_g_gene_level_var_prior[k]

        self.alpha_g_phi_hyp_prior = alpha_g_phi_hyp_prior
        self.w_sf_mean_var_ratio = w_sf_mean_var_ratio
        self.gene_add_alpha_hyp_prior = gene_add_alpha_hyp_prior
        self.gene_add_mean_hyp_prior = gene_add_mean_hyp_prior

        factors_per_groups = A_factors_per_location / Y_groups_per_location

        # compute hyperparameters from mean and sd
        self.m_g_gene_level_prior = m_g_gene_level_prior
        self.register_buffer(
            "m_g_shape",
            torch.tensor(
                (self.m_g_gene_level_prior["mean"] ** 2)
                / (self.m_g_gene_level_prior["sd"] ** 2)
            ),
        )
        self.register_buffer(
            "m_g_rate",
            torch.tensor(
                self.m_g_gene_level_prior["mean"]
                / (self.m_g_gene_level_prior["sd"] ** 2)
            ),
        )
        self.register_buffer(
            "m_g_mean_var", torch.tensor(self.m_g_gene_level_prior["mean_var_ratio"])
        )

        self.cell_state_mat = cell_state_mat
        self.register_buffer("cell_state", torch.tensor(cell_state_mat.T))

        self.register_buffer("N_cells_per_location", torch.tensor(N_cells_per_location))
        self.register_buffer("factors_per_groups", torch.tensor(factors_per_groups))
        self.register_buffer(
            "Y_groups_per_location", torch.tensor(Y_groups_per_location)
        )
        self.register_buffer(
            "N_cells_mean_var_ratio", torch.tensor(N_cells_mean_var_ratio)
        )

        self.register_buffer(
            "alpha_g_phi_hyp_prior_alpha",
            torch.tensor(self.alpha_g_phi_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "alpha_g_phi_hyp_prior_beta",
            torch.tensor(self.alpha_g_phi_hyp_prior["beta"]),
        )
        self.register_buffer(
            "gene_add_alpha_hyp_prior_alpha",
            torch.tensor(self.gene_add_alpha_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "gene_add_alpha_hyp_prior_beta",
            torch.tensor(self.gene_add_alpha_hyp_prior["beta"]),
        )
        self.register_buffer(
            "gene_add_mean_hyp_prior_alpha",
            torch.tensor(self.gene_add_mean_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "gene_add_mean_hyp_prior_beta",
            torch.tensor(self.gene_add_mean_hyp_prior["beta"]),
        )

        self.register_buffer(
            "w_sf_mean_var_ratio_tensor", torch.tensor(self.w_sf_mean_var_ratio)
        )

        self.register_buffer("n_factors_tensor", torch.tensor(self.n_factors))
        self.register_buffer("n_groups_tensor", torch.tensor(self.n_groups))

        self.register_buffer("ones", torch.ones((1, 1)))
        self.register_buffer("ones_1_n_groups", torch.ones((1, self.n_groups)))
        self.register_buffer("eps", torch.tensor(1e-8))

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x_data = tensor_dict[_CONSTANTS.X_KEY]
        ind_x = tensor_dict["ind_x"].long().squeeze()
        batch_index = tensor_dict[_CONSTANTS.BATCH_KEY]
        return (x_data, ind_x, batch_index), {}

    def _get_fn_args_full_data(self, adata):
        x_data = get_from_registry(adata, _CONSTANTS.X_KEY)
        if issparse(x_data):
            x_data = np.asarray(x_data.toarray())
        x_data = torch.tensor(x_data.astype("float32"))
        ind_x = torch.tensor(get_from_registry(adata, "ind_x"))
        batch_index = torch.tensor(get_from_registry(adata, _CONSTANTS.BATCH_KEY))
        return (x_data, ind_x, batch_index), {}

    def create_plates(self, x_data, idx, batch_index):
        return pyro.plate("obs_plate", size=self.n_obs, dim=-2, subsample=idx)

    def list_obs_plate_vars(self):
        """Create a dictionary with the name of observation/minibatch plate,
        indexes of model args to provide to encoder,
        variable names that belong to the observation plate
        and the number of dimensions in non-plate axis of each variable"""

        return {
            "name": "obs_plate",
            "in": [0],  # expression data + (optional) batch index
            "sites": {
                "n_s_cells_per_location": 1,
                "y_s_groups_per_location": 1,
                "z_sr_groups_factors": self.n_groups,
                "w_sf": self.n_factors,
                "l_s_add": 1,
            },
        }

    def forward(self, x_data, idx, batch_index):

        obs2sample = one_hot(batch_index, self.n_batch)

        obs_plate = self.create_plates(x_data, idx, batch_index)

        # =====================Gene expression level scaling m_g======================= #
        # Explains difference in sensitivity for each gene between single cell and spatial technology

        m_g_alpha_hyp = pyro.sample(
            "m_g_alpha_hyp",
            dist.Gamma(self.m_g_shape * self.m_g_mean_var, self.m_g_mean_var),
        )

        m_g_beta_hyp = pyro.sample(
            "m_g_beta_hyp",
            dist.Gamma(self.m_g_rate * self.m_g_mean_var, self.m_g_mean_var),
        )

        m_g = pyro.sample(
            "m_g",
            dist.Gamma(m_g_alpha_hyp, m_g_beta_hyp)
            .expand([1, self.n_vars])
            .to_event(2),
        )

        # =====================Cell abundances w_sf======================= #
        # factorisation prior on w_sf models similarity in locations
        # across cell types f and reflects the absolute scale of w_sf
        with obs_plate:
            n_s_cells_per_location = pyro.sample(
                "n_s_cells_per_location",
                dist.Gamma(
                    self.N_cells_per_location * self.N_cells_mean_var_ratio,
                    self.N_cells_mean_var_ratio,
                ),
            )

            y_s_groups_per_location = pyro.sample(
                "y_s_groups_per_location",
                dist.Gamma(self.Y_groups_per_location, self.ones),
            )

        # cell group loadings
        shape = self.ones_1_n_groups * y_s_groups_per_location / self.n_groups_tensor
        rate = self.ones_1_n_groups / (n_s_cells_per_location / y_s_groups_per_location)
        with obs_plate:
            z_sr_groups_factors = pyro.sample(
                "z_sr_groups_factors",
                dist.Gamma(
                    shape, rate
                ),  # .to_event(1)#.expand([self.n_groups]).to_event(1)
            )  # (n_obs, n_groups)

        k_r_factors_per_groups = pyro.sample(
            "k_r_factors_per_groups",
            dist.Gamma(self.factors_per_groups, self.ones)
            .expand([self.n_groups, 1])
            .to_event(2),
        )  # (self.n_groups, 1)

        c2f_shape = k_r_factors_per_groups / self.n_factors_tensor

        x_fr_group2fact = pyro.sample(
            "x_fr_group2fact",
            dist.Gamma(c2f_shape, k_r_factors_per_groups)
            .expand([self.n_groups, self.n_factors])
            .to_event(2),
        )  # (self.n_groups, self.n_factors)

        with obs_plate:
            w_sf_mu = z_sr_groups_factors @ x_fr_group2fact
            w_sf = pyro.sample(
                "w_sf",
                dist.Gamma(
                    w_sf_mu * self.w_sf_mean_var_ratio_tensor,
                    self.w_sf_mean_var_ratio_tensor,
                ),
            )  # (self.n_obs, self.n_factors)

        # =====================Location-specific additive component======================= #
        l_s_add_alpha = pyro.sample("l_s_add_alpha", dist.Gamma(self.ones, self.ones))
        l_s_add_beta = pyro.sample("l_s_add_beta", dist.Gamma(self.ones, self.ones))

        with obs_plate:
            l_s_add = pyro.sample(
                "l_s_add", dist.Gamma(l_s_add_alpha, l_s_add_beta)
            )  # (self.n_obs, 1)

        # =====================Gene-specific additive component ======================= #
        # per gene molecule contribution that cannot be explained by
        # cell state signatures (e.g. background, free-floating RNA)
        s_g_gene_add_alpha_hyp = pyro.sample(
            "s_g_gene_add_alpha_hyp",
            dist.Gamma(
                self.gene_add_alpha_hyp_prior_alpha, self.gene_add_alpha_hyp_prior_beta
            ),
        )
        s_g_gene_add_mean = pyro.sample(
            "s_g_gene_add_mean",
            dist.Gamma(
                self.gene_add_mean_hyp_prior_alpha,
                self.gene_add_mean_hyp_prior_beta,
            )
            .expand([self.n_batch, 1])
            .to_event(2),
        )  # (self.n_batch)
        s_g_gene_add_alpha_e_inv = pyro.sample(
            "s_g_gene_add_alpha_e_inv",
            dist.Exponential(s_g_gene_add_alpha_hyp)
            .expand([self.n_batch, 1])
            .to_event(2),
        )  # (self.n_batch)
        s_g_gene_add_alpha_e = self.ones / s_g_gene_add_alpha_e_inv.pow(2)

        s_g_gene_add = pyro.sample(
            "s_g_gene_add",
            dist.Gamma(s_g_gene_add_alpha_e, s_g_gene_add_alpha_e / s_g_gene_add_mean)
            .expand([self.n_batch, self.n_vars])
            .to_event(2),
        )  # (self.n_batch, n_vars)

        # =====================Gene-specific overdispersion ======================= #
        alpha_g_phi_hyp = pyro.sample(
            "alpha_g_phi_hyp",
            dist.Gamma(
                self.alpha_g_phi_hyp_prior_alpha, self.alpha_g_phi_hyp_prior_beta
            ),
        )
        alpha_g_inverse = pyro.sample(
            "alpha_g_inverse",
            dist.Exponential(alpha_g_phi_hyp)
            .expand([self.n_batch, self.n_vars])
            .to_event(2),
        )  # (self.n_batch, self.n_vars)

        # =====================Expected expression ======================= #
        # expected expression
        mu = (w_sf @ self.cell_state) * m_g + (obs2sample @ s_g_gene_add) + l_s_add
        theta = obs2sample @ (self.ones / alpha_g_inverse.pow(2))
        # convert mean and overdispersion to total count and logits
        total_count, logits = _convert_mean_disp_to_counts_logits(
            mu, theta, eps=self.eps
        )

        # =====================DATA likelihood ======================= #
        # Likelihood (sampling distribution) of data_target & add overdispersion via NegativeBinomial
        with obs_plate:
            pyro.sample(
                "data_target",
                dist.NegativeBinomial(total_count=total_count, logits=logits),
                # NegativeBinomial(mu=mu, theta=theta),
                obs=x_data,
            )

        # =====================Compute mRNA count from each factor in locations  ======================= #
        mRNA = w_sf * (self.cell_state * m_g).sum(-1)
        pyro.deterministic("u_sf_mRNA_factors", mRNA)

    def compute_expected(self, samples, adata):
        r"""Compute expected expression of each gene in each location. Useful for evaluating how well
        the model learned expression pattern of all genes in the data.
        """
        obs2sample = get_from_registry(adata, _CONSTANTS.BATCH_KEY)
        obs2sample = pd.get_dummies(obs2sample.flatten())
        mu = (
            np.dot(samples["w_sf"], self.cell_state_mat.T) * samples["m_g"]
            + np.dot(obs2sample, samples["s_g_gene_add"])
            + samples["l_s_add"]
        )
        alpha = np.dot(obs2sample, 1 / np.power(samples["alpha_g_inverse"], 2))

        return {"mu": mu, "alpha": alpha}
