import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, Beta, Gamma, kl_divergence as kl


from scvi.models.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.models.modules import Encoder, DecoderSCVI, LinearDecoderSCVI
from scvi.models.vae import VAE
from scipy.special import logit
from scvi.models.utils import one_hot

torch.backends.cudnn.benchmark = True

class AutoZIVAE(VAE):

    def __init__(
        self,
        n_input: int,
        alpha_prior: float = 0.5,
        beta_prior: float = 0.5,
        minimal_dropout: float = 0.01,
        **args,
    ):

        if 'reconstruction_loss' in args:
            raise ValueError(
                "No reconstruction loss must be specified for AutoZI : it is 'autozinb'."
            )

        super().__init__(n_input, **args)
        self.reconstruction_loss = 'autozinb'
        self.minimal_dropout = minimal_dropout


        # Parameters of prior Bernoulli Beta distribution : alpha + beta = 1 if only one is specified
        if beta_prior is None and alpha_prior is not None:
            beta_prior = 1. - alpha_prior
        if alpha_prior is None and beta_prior is not None:
            alpha_prior = 1. - beta_prior

        device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

        # Create parameters for Bernoulli Beta prior and posterior distributions
        # Each paramer, whose values are in (0,1), is encoded as its logit, in the set of real numbers

        if self.dispersion == "gene":
            self.alpha_posterior_logit = torch.nn.Parameter(torch.randn(n_input, ))
            self.beta_posterior_logit = torch.nn.Parameter(torch.randn(n_input, ))
            self.alpha_prior_logit = torch.nn.Parameter(torch.randn(1, )) \
                if alpha_prior is None \
                else torch.Tensor([logit(alpha_prior)]).to(device)
            self.beta_prior_logit = torch.nn.Parameter(torch.randn(1, )) \
                if beta_prior is None \
                else torch.Tensor([logit(beta_prior)]).to(device)

        elif self.dispersion == "gene-batch":
            self.alpha_posterior_logit = torch.nn.Parameter(torch.randn(n_input, self.n_batch))
            self.beta_posterior_logit = torch.nn.Parameter(torch.randn(n_input, self.n_batch))
            self.alpha_prior_logit = torch.nn.Parameter(torch.randn(1, self.n_batch)) \
                if alpha_prior is None \
                else torch.Tensor([logit(alpha_prior)]).to(device)
            self.beta_prior_logit = torch.nn.Parameter(torch.randn(1, self.n_batch)) \
                if beta_prior is None \
                else torch.Tensor([logit(beta_prior)]).to(device)

        elif self.dispersion == "gene-label":
            self.alpha_posterior_logit = torch.nn.Parameter(torch.randn(n_input, self.n_labels))
            self.beta_posterior_logit = torch.nn.Parameter(torch.randn(n_input, self.n_labels))
            self.alpha_prior_logit = torch.nn.Parameter(torch.randn(1, self.n_labels)) \
                if alpha_prior is None \
                else torch.Tensor([logit(alpha_prior)]).to(device)
            self.beta_prior_logit = torch.nn.Parameter(torch.randn(1, self.n_labels)) \
                if beta_prior is None \
                else torch.Tensor([logit(beta_prior)]).to(device)

        else:  # gene-cell
            raise Exception("Gene-cell not implemented yet for AutoZI")




    def get_alphas_betas(
        self,
        as_numpy=True,
    ):
        # Return parameters of Bernoulli Beta distributions in a dictionary

        outputs = {}
        outputs['alpha_posterior'] = torch.sigmoid(self.alpha_posterior_logit)
        outputs['beta_posterior'] = torch.sigmoid(self.beta_posterior_logit)
        outputs['alpha_prior'] = torch.sigmoid(self.alpha_prior_logit)
        outputs['beta_prior'] = torch.sigmoid(self.beta_prior_logit)

        if as_numpy:
            for key,value in outputs.items():
                outputs[key] = value.detach().cpu().numpy() if value.requires_grad else value.cpu().numpy()

        return outputs



    def sample_from_beta_distribution(self, alpha, beta, eps_gamma=1e-30, eps_sample=1e-7):
        # Sample from a Beta distribution using the reparameterization trick.
        # Problem : it is not implemented in CUDA yet
        # Workaround : sample X and Y from Gamma(alpha,1) and Gamma(beta,1), the Beta sample is X/(X+Y)
        # Warning : use logs and perform logsumexp to avoid numerical issues

        # Sample from Gamma
        sample_x_log = torch.log(Gamma(alpha, 1).rsample() + eps_gamma)
        sample_y_log = torch.log(Gamma(beta, 1).rsample() + eps_gamma)

        # Sum using logsumexp (note : eps_gamma is used to prevent numerical issues with perfect
        # 0 and 1 final Beta samples
        sample_xy_log_max = torch.max(sample_x_log, sample_y_log)
        sample_xplusy_log = sample_xy_log_max \
                                       + torch.log(torch.exp(sample_x_log - sample_xy_log_max) \
                                                   + torch.exp(sample_y_log - sample_xy_log_max))
        sample_log = sample_x_log - sample_xplusy_log
        sample = eps_sample + (1 - 2*eps_sample)*torch.exp(sample_log)

        return sample


    def get_reconstruction_loss(self, x, px_rate, px_r, px_dropout, bernoulli_params, eps_log=1e-8):

        # LLs for NB and ZINB
        ll_zinb = torch.log(1. - bernoulli_params + eps_log) + log_zinb_positive(x, px_rate, px_r, px_dropout,\
                                                                        return_gene_specific=True)
        ll_nb = torch.log(bernoulli_params + eps_log) + log_nb_positive(x, px_rate, px_r, return_gene_specific=True)

        # Reconstruction loss using a logsumexp-type computation
        ll_max = torch.max(ll_zinb, ll_nb)
        ll_tot = ll_max + torch.log(torch.exp(ll_nb-ll_max) + torch.exp(ll_zinb-ll_max))
        reconst_loss = -ll_tot.sum(dim=-1)

        return reconst_loss


    def rescale_bernoulli_dispersion(self, bernoulli_params, x, batch_index=None, y=None):
        if self.dispersion == "gene-label":
            one_hot_label = one_hot(y, self.n_labels)
            # If we sampled several random Bernoulli parameters
            if len(bernoulli_params.shape) == 2:
                bernoulli_params = F.linear(one_hot_label, bernoulli_params)
            else:
                bernoulli_params_res = []
                for sample in range(bernoulli_params.shape[0]):
                    bernoulli_params_res.append(F.linear(one_hot_label, bernoulli_params[sample]))
                bernoulli_params = torch.stack(bernoulli_params_res)
        elif self.dispersion == "gene-batch":
            one_hot_batch = one_hot(batch_index, self.n_batch)
            if len(bernoulli_params.shape) == 2:
                bernoulli_params = F.linear(one_hot_batch, bernoulli_params)
            # If we sampled several random Bernoulli parameters
            else:
                bernoulli_params_res = []
                for sample in range(bernoulli_params.shape[0]):
                    bernoulli_params_res.append(F.linear(one_hot_batch, bernoulli_params[sample]))
                bernoulli_params = torch.stack(bernoulli_params_res)

        return bernoulli_params



    def inference(self, x, batch_index=None, y=None, n_samples=1, eps_log=1e-8):
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        # Sampling
        qz_m, qz_v, z = self.z_encoder(x_, y)
        ql_m, ql_v, library = self.l_encoder(x_)


        outputs = self.get_alphas_betas(as_numpy=False)
        alpha_posterior = outputs['alpha_posterior']
        beta_posterior = outputs['beta_posterior']
        alpha_prior = outputs['alpha_prior']
        beta_prior = outputs['beta_prior']

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            z = Normal(qz_m, qz_v.sqrt()).sample()
            ql_m = ql_m.unsqueeze(0).expand((n_samples, ql_m.size(0), ql_m.size(1)))
            ql_v = ql_v.unsqueeze(0).expand((n_samples, ql_v.size(0), ql_v.size(1)))
            library = Normal(ql_m, ql_v.sqrt()).sample()
            alpha_posterior = alpha_posterior.unsqueeze(0).expand((n_samples, alpha_posterior.size(0))) \
                if self.dispersion == 'gene' else \
                alpha_posterior.unsqueeze(0).expand((n_samples, alpha_posterior.size(0), alpha_posterior.size(1)))
            beta_posterior = beta_posterior.unsqueeze(0).expand((n_samples, beta_posterior.size(0))) \
                if self.dispersion == 'gene' else \
                beta_posterior.unsqueeze(0).expand((n_samples, beta_posterior.size(0), beta_posterior.size(1)))

        bernoulli_params = self.sample_from_beta_distribution(alpha_posterior, beta_posterior)
        bernoulli_params = self.rescale_bernoulli_dispersion(bernoulli_params, x, batch_index, y)

        px_scale, px_r, px_rate, px_dropout = self.decoder(
            self.dispersion, z, library, batch_index, y
        )
        if self.dispersion == "gene-label":
            px_r = F.linear(
                one_hot(y, self.n_labels), self.px_r
            )  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r
        px_r = torch.exp(px_r)

        # Rescale dropout
        if self.minimal_dropout > 0.:
            dropout_prob_rescaled = self.minimal_dropout + (1. - self.minimal_dropout) * torch.sigmoid(px_dropout)
            px_dropout_rescaled = torch.log(dropout_prob_rescaled / (1. - dropout_prob_rescaled + eps_log))
        else:
            px_dropout_rescaled = px_dropout

        return dict(
            px_scale=px_scale,
            px_r=px_r,
            px_rate=px_rate,
            px_dropout=px_dropout_rescaled,
            qz_m=qz_m,
            qz_v=qz_v,
            z=z,
            ql_m=ql_m,
            ql_v=ql_v,
            library=library,
            alpha_posterior=alpha_posterior,
            beta_posterior=beta_posterior,
            alpha_prior=alpha_prior,
            beta_prior=beta_prior,
            bernoulli_params=bernoulli_params,
        )

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        r""" Returns the reconstruction loss and the Kullback divergences

        :param x: tensor of values with shape (batch_size, n_input)
        :param local_l_mean: tensor of means of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param local_l_var: tensor of variancess of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape (batch_size, n_labels)
        :return: the reconstruction loss and the Kullback divergences
        :rtype: 2-tuple of :py:class:`torch.FloatTensor`
        """
        # Parameters for z latent distribution
        outputs = self.inference(x, batch_index, y)
        qz_m = outputs["qz_m"]
        qz_v = outputs["qz_v"]
        ql_m = outputs["ql_m"]
        ql_v = outputs["ql_v"]
        px_rate = outputs["px_rate"]
        px_r = outputs["px_r"]
        px_dropout = outputs["px_dropout"]
        alpha_posterior = outputs["alpha_posterior"]
        beta_posterior = outputs["beta_posterior"]
        alpha_prior = outputs["alpha_prior"]
        beta_prior = outputs["beta_prior"]
        bernoulli_params = outputs["bernoulli_params"]

        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(
            dim=1
        )
        kl_divergence_l = kl(
            Normal(ql_m, torch.sqrt(ql_v)),
            Normal(local_l_mean, torch.sqrt(local_l_var)),
        ).sum(dim=1)

        kl_divergence_bernoulli = kl(Beta(alpha_posterior, beta_posterior), Beta(alpha_prior, beta_prior)).sum()

        reconst_loss = self.get_reconstruction_loss(x, px_rate, px_r, px_dropout, bernoulli_params)

        return reconst_loss + kl_divergence_l, kl_divergence_z, kl_divergence_bernoulli


