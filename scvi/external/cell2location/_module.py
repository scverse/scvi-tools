import os

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
from pyro.infer.autoguide import AutoNormal, init_to_mean
from pyro.nn import PyroModule, pyro_method

from scvi import _CONSTANTS
from scvi.compose import PyroBaseModuleClass
from scvi.distributions._negative_binomial import _convert_mean_disp_to_counts_logits

from torch.distributions import constraints
from torch.distributions.transforms import Transform
from torch.distributions.constraint_registry import biject_to, transform_to

class SoftplusTransform(Transform):
    r"""
    Transform via the mapping :math:`\text{Softplus}(x) = \log(1 + \exp(x))`.
    """
    domain = constraints.real_vector
    codomain = constraints.positive
    
    def __eq__(self, other):
        return isinstance(other, SoftplusTransform)

    def _call(self, x):
        return (1+x.exp()).log()

    def _inverse(self, y):
        return y.expm1().log()
      
    def log_abs_det_jacobian(self, x, y):
        return -(1+(-x).exp()).log()

# replace exp transform with softplus https://github.com/pyro-ppl/numpyro/issues/855
class _Positive(constraints.Constraint):
    
    def __init__(self):
        
        self.lower_bound = 0
        super().__init__()
    
    def check(self, value):
        return self.lower_bound < value

    def __repr__(self):
        fmt_string = self.__class__.__name__[1:]
        fmt_string += '(lower_bound={})'.format(self.lower_bound)
        return fmt_string

@biject_to.register(_Positive)
@transform_to.register(_Positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()

class exp_to_softplus(pyro.poutine.messenger.Messenger):
    
    def process_message(self, msg):
        if msg["type"] == "param" and msg["name"].endswith("_scale"):
            msg["kwargs"]["constraint"] = _Positive()

# Define helper gamma distribution
def Gamma(mu=None, sigma=None, alpha=None, beta=None, shape=None):
    r"""
    Function that converts mu/sigma Gamma distribution parametrisation into alpha/beta and
    returns pyro Gamma distribution with an event of a given shape.
    A thin wrapper over pyro.dist.Gamma.

    :param mu: mean of Gamma distribution
    :param sigma: variance of Gamma distribution. Note!!! this is :math:`\sigma^2` not :math:`\sigma`
    :param alpha: shape parameter of Gamma distribution (mu ** 2 / sigma)
    :param beta: rate parameter of Gamma distribution (mu / sigma)
    :param shape: shape of the event / resulting variable. When None the shape is guessed based on input parameters

    :return: pyro Gamma distribution class
    """
    if alpha is not None and beta is not None:
        pass
    elif mu is not None and sigma is not None:
        alpha = mu ** 2 / sigma
        beta = mu / sigma
    else:
        raise ValueError('Define (mu and sigma_sq) or (alpha and beta).')

    if shape is None:
        alpha = torch.tensor(alpha)
        beta = torch.tensor(beta)
    else:
        alpha = torch.ones(shape) * torch.tensor(alpha)
        beta = torch.ones(shape) * torch.tensor(beta)
    return dist.Gamma(alpha, beta)

class LocationModelLinearDependentWMultiExperiment(PyroBaseModuleClass):
    
    def __init__(self, n_obs, n_var, n_fact, n_exper, batch_size,
                 cell_state_mat, n_comb: int = 50,
                 m_g_gene_level_prior={'mean': 1 / 2, 'sd': 1 / 4},
                 m_g_gene_level_var_prior={'mean_var_ratio': 1},
                 cell_number_prior={'N_cells_per_location': 8,
                               'A_factors_per_location': 7,
                               'Y_combs_per_location': 2.5},
                 cell_number_var_prior={'N_cells_mean_var_ratio': 1,
                                   'A_factors_mean_var_ratio': 1,
                                   'Y_combs_mean_var_ratio': 1},
                 alpha_g_phi_hyp_prior={'mean': 3, 'sd': 1},
                 gene_add_alpha_hyp_prior={'mean': 3, 'sd': 1},
                 gene_add_mean_hyp_prior={'alpha': 1, 'beta': 100},
                 w_sf_mean_var_ratio=5):

        super().__init__()
        
        self.n_obs = n_obs
        self.n_var = n_var
        self.n_fact = n_fact
        self.n_exper = n_exper
        self.batch_size = batch_size
        self.n_comb = n_comb
        
        for k in m_g_gene_level_var_prior.keys():
            m_g_gene_level_prior[k] = m_g_gene_level_var_prior[k]
        for k in m_g_gene_level_prior.keys():
            m_g_gene_level_prior[k] = np.array(m_g_gene_level_prior[k]).reshape((1, 1))

        self.m_g_gene_level_prior = m_g_gene_level_prior
        self.alpha_g_phi_hyp_prior = alpha_g_phi_hyp_prior
        self.w_sf_mean_var_ratio = w_sf_mean_var_ratio
        self.gene_add_alpha_hyp_prior = gene_add_alpha_hyp_prior
        self.gene_add_mean_hyp_prior = gene_add_mean_hyp_prior
        
        cell_number_prior['factors_per_combs'] = (cell_number_prior['A_factors_per_location'] /
                                                  cell_number_prior['Y_combs_per_location'])
        for k in cell_number_var_prior.keys():
            cell_number_prior[k] = cell_number_var_prior[k]
        for k in cell_number_prior.keys():
            cell_number_prior[k] = np.array(cell_number_prior[k]).reshape((1, 1))
        self.cell_number_prior = cell_number_prior
        
        self.cell_state_mat = cell_state_mat
        self.register_buffer("cell_state", torch.tensor(cell_state_mat.T))
        
        self.guide = AutoNormal(self.model, init_loc_fn=init_to_mean, 
                                create_plates=self.create_plates)
        # replace exp transform with softplus https://github.com/pyro-ppl/numpyro/issues/855
        #self.guide = exp_to_softplus(self.guide)

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x_data = tensor_dict[_CONSTANTS.X_KEY]
        tensor_dict["ind_x"].squeeze().type(torch.cuda.LongTensor)
        obs2sample = tensor_dict["obs2sample"]
        return (x_data, ind_x, obs2sample), {}
    
    
    def create_plates(self, x_data, idx, obs2sample):
        return [pyro.plate("obs_axis", self.n_obs, dim=-2, 
                           subsample_size=self.batch_size, 
                           subsample=idx),
                pyro.plate("var_axis", self.n_var, dim=-1),
                pyro.plate("factor_axis", self.n_fact, dim=-1),
                pyro.plate("combination_axis", self.n_comb, dim=-3),
                pyro.plate("experim_axis", self.n_exper, dim=-2)]

    def model(self, x_data, idx, obs2sample):
        # register module with Pyro
        pyro.module("cell2location", self)
        
        obs_axis, var_axis, factor_axis, combination_axis, experim_axis = self.create_plates(x_data, idx, obs2sample)
        
        # =====================Gene expression level scaling m_g======================= #
        # Explains difference in sensitivity for each gene between single cell and spatial technology
        # compute hyperparameters from mean and sd
        shape = self.m_g_gene_level_prior['mean'] ** 2 / self.m_g_gene_level_prior['sd'] ** 2
        rate = self.m_g_gene_level_prior['mean'] / self.m_g_gene_level_prior['sd'] ** 2
        shape_var = shape / self.m_g_gene_level_prior['mean_var_ratio']
        rate_var = rate / self.m_g_gene_level_prior['mean_var_ratio']
        
        m_g_alpha_hyp = pyro.sample('m_g_alpha_hyp',
                                                Gamma(mu=shape,
                                                      sigma=shape_var,
                                                      shape=(1, 1)))

        m_g_beta_hyp = pyro.sample('m_g_beta_hyp',
                                               Gamma(mu=rate,
                                                     sigma=rate_var,
                                                     shape=(1, 1)))
        with var_axis:
            m_g = pyro.sample('m_g', Gamma(alpha=m_g_alpha_hyp,
                                                beta=m_g_beta_hyp,
                                                shape=(1, self.n_var)))

        # =====================Cell abundances w_sf======================= #
        # factorisation prior on w_sf models similarity in locations 
        # between cell types f and reflects the absolute scale of w_sf
        with obs_axis:
            N_s_cells_per_location = pyro.sample('N_s_cells_per_location',
                                              Gamma(mu=self.cell_number_prior['N_cells_per_location'],
                                                    sigma=self.cell_number_prior['N_cells_per_location'] \
                                                                  / self.cell_number_prior['N_cells_mean_var_ratio'],
                                                    shape=None))

            Y_s_combs_per_location = pyro.sample('Y_s_combs_per_location',
                                             Gamma(mu=self.cell_number_prior['Y_combs_per_location'],
                                                   sigma=self.cell_number_prior['Y_combs_per_location'] \
                                                                 / self.cell_number_prior['Y_combs_mean_var_ratio'],
                                                   shape=None))
            
            with combination_axis:
                shape = Y_s_combs_per_location / self.n_comb
                rate = torch.ones([self.n_comb, 1, 1]) / N_s_cells_per_location * Y_s_combs_per_location
                z_sr_combs_factors = pyro.sample('z_sr_combs_factors',
                                                 Gamma(alpha=shape,
                                                       beta=rate,
                                                       shape=None)) # (n_comb, n_obs, 1)
        with combination_axis:
            K_r_factors_per_combs = pyro.sample('K_r_factors_per_combs',
                                                 Gamma(mu=self.cell_number_prior['factors_per_combs'],
                                                       sigma=self.cell_number_prior['factors_per_combs'] \
                                                                     / self.cell_number_prior['A_factors_mean_var_ratio'],
                                                       shape=(self.n_comb, 1, 1)))

            c2f_shape = K_r_factors_per_combs / self.n_fact
            
        with factor_axis, combination_axis:
            x_fr_comb2fact = pyro.sample('x_fr_comb2fact',
                                             Gamma(alpha=c2f_shape,
                                                   beta=K_r_factors_per_combs,
                                                   shape=(self.n_comb, 1, self.n_fact)))

        with obs_axis:
            with factor_axis:
                w_sf_mu = z_sr_combs_factors.squeeze(-1).T @ x_fr_comb2fact.squeeze(-2)
                w_sf_sigma = w_sf_mu / self.w_sf_mean_var_ratio
                w_sf = pyro.sample('w_sf', Gamma(mu=w_sf_mu, sigma=w_sf_sigma)) # (self.n_obs, self.n_fact)

        # =====================Location-specific additive component======================= #
        l_s_add_alpha = pyro.sample('l_s_add_alpha', Gamma(alpha=1, beta=1, shape=(1, 1)))
        l_s_add_beta = pyro.sample('l_s_add_beta', Gamma(alpha=1, beta=1, shape=(1, 1)))
                                        
        with obs_axis:
            l_s_add = pyro.sample('l_s_add', Gamma(alpha=l_s_add_alpha,
                                                        beta=l_s_add_beta,
                                                        shape=None)) # (self.n_obs, 1)

        # =====================Gene-specific additive component ======================= #
        # per gene molecule contribution that cannot be explained by 
        # cell state signatures (e.g. background, free-floating RNA)
        with experim_axis:
            s_g_gene_add_mean = pyro.sample('s_g_gene_add_mean',
                                      Gamma(alpha=self.gene_add_mean_hyp_prior['alpha'],
                                      beta=self.gene_add_mean_hyp_prior['beta'], shape=(self.n_exper, 1)))
            s_g_gene_add_alpha_hyp = pyro.sample('s_g_gene_add_alpha_hyp',
                                           Gamma(mu=self.gene_add_alpha_hyp_prior['mean'],
                                           sigma=self.gene_add_alpha_hyp_prior['sd'], shape=(1, 1)))
            s_g_gene_add_alpha_e_inv = pyro.sample('s_g_gene_add_alpha_e_inv',
                                                   dist.Exponential(s_g_gene_add_alpha_hyp * torch.ones([self.n_exper, 1])))
            s_g_gene_add_alpha_e = torch.ones((1, 1)) / torch.pow(s_g_gene_add_alpha_e_inv, 2)
            with var_axis:
                s_g_gene_add = pyro.sample('s_g_gene_add', Gamma(s_g_gene_add_alpha_e,
                                                   s_g_gene_add_alpha_e / s_g_gene_add_mean,
                                                  shape=(self.n_exper, self.n_var)))

        # =====================Gene-specific overdispersion ======================= #
        alpha_g_phi_hyp = pyro.sample('alpha_g_phi_hyp',
                                   Gamma(mu=self.alpha_g_phi_hyp_prior['mean'],
                                         sigma=self.alpha_g_phi_hyp_prior['sd'],
                                         shape=(1, 1)))
        with experim_axis:
            with var_axis:
                alpha_g_inverse = pyro.sample('alpha_g_inverse', 
                                              dist.Exponential(alpha_g_phi_hyp 
                                                               * torch.ones([self.n_exper, self.n_var])))

        # =====================Expected expression ======================= #
        # expected expression
        #print(w_sf_mu.shape)
        #print(z_sr_combs_factors.squeeze(-1).T.shape)
        #print(x_fr_comb2fact.squeeze(-2).shape)
        #print(w_sf.shape)
        #print(self.cell_state.shape)
        #print(m_g.shape)
        #print(obs2sample.shape)
        #print(s_g_gene_add.shape)
        #print(l_s_add.shape)
        #print(alpha_g_inverse.shape)
        mu = torch.matmul(w_sf, self.cell_state) * m_g + torch.matmul(obs2sample, s_g_gene_add) + l_s_add
        theta = torch.matmul(obs2sample, torch.ones([1, 1]) / (alpha_g_inverse * alpha_g_inverse))
        # convert mean and overdispersion to total count and logits (input to NB from pyro)
        total_count, logits = _convert_mean_disp_to_counts_logits(mu, theta, eps=1e-8)

        # =====================DATA likelihood ======================= #
        # Likelihood (sampling distribution) of data_target & add overdispersion via NegativeBinomial
        with var_axis:
            with obs_axis:
                self.data_target = pyro.sample('data_target',
                                               dist.NegativeBinomial(total_count=total_count,
                                                                     logits=logits),
                                               obs=x_data)

        # =====================Compute mRNA count from each factor in locations  ======================= #
        print(idx)
        mRNA = (w_sf * (self.cell_state * m_g).sum(-1))
        u_sf_mRNA_factors = pyro.deterministic('u_sf_mRNA_factors', mRNA)

    def compute_expected(self, obs2sample):
        r"""Compute expected expression of each gene in each location. Useful for evaluating how well
            the model learned expression pattern of all genes in the data.
        """

        # compute the poisson rate
        self.mu = (np.dot(self.samples['post_sample_means']['w_sf'],
                          self.cell_state_mat.T)
                   * self.samples['post_sample_means']['m_g']
                   + np.dot(obs2sample,
                            self.samples['post_sample_means']['s_g_gene_add'])
                   + self.samples['post_sample_means']['l_s_add'])
        self.alpha = np.dot(obs2sample,
                            1 / (self.samples['post_sample_means']['alpha_g_inverse'] * self.samples['post_sample_means']['alpha_g_inverse']))

     
 # set default module
Cell2locationModule = LocationModelLinearDependentWMultiExperiment
