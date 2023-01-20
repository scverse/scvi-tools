from typing import Dict, Literal, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn.functional as F
from scipy.special import logit
from torch.distributions import Beta, Gamma, Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.module.base import LossRecorder, auto_move_data
from scvi.nn import one_hot

from ._vae import VAE

torch.backends.cudnn.benchmark = True


class AutoZIVAE(VAE):
    """
    Implementation of the AutoZI model :cite:p:`Clivio19`.

    Parameters
    ----------
    n_input
        Number of input genes
    alpha_prior
        Float denoting the alpha parameter of the prior Beta distribution of
        the zero-inflation Bernoulli parameter. Should be between 0 and 1, not included.
        When set to ``None``, will be set to 1 - beta_prior if beta_prior is not ``None``,
        otherwise the prior Beta distribution will be learned on an Empirical Bayes fashion.
    beta_prior
        Float denoting the beta parameter of the prior Beta distribution of
        the zero-inflation Bernoulli parameter. Should be between 0 and 1, not included.
        When set to ``None``, will be set to 1 - alpha_prior if alpha_prior is not ``None``,
        otherwise the prior Beta distribution will be learned on an Empirical Bayes fashion.
    minimal_dropout
        Float denoting the lower bound of the cell-gene ZI rate in the ZINB component.
        Must be non-negative. Can be set to 0 but not recommended as this may make
        the mixture problem ill-defined.
    zero_inflation: One of the following

        * ``'gene'`` - zero-inflation Bernoulli parameter of AutoZI is constant per gene across cells
        * ``'gene-batch'`` - zero-inflation Bernoulli parameter can differ between different batches
        * ``'gene-label'`` - zero-inflation Bernoulli parameter can differ between different labels
        * ``'gene-cell'`` - zero-inflation Bernoulli parameter can differ for every gene in every cell


    See VAE docstring (scvi/models/vae.py) for more parameters. ``reconstruction_loss`` should not be specified.

    Examples
    --------
    >>> gene_dataset = CortexDataset()
    >>> autozivae = AutoZIVAE(gene_dataset.nb_genes, alpha_prior=0.5, beta_prior=0.5, minimal_dropout=0.01)

    """

    def __init__(
        self,
        n_input: int,
        alpha_prior: Tunable[float] = 0.5,
        beta_prior: Tunable[float] = 0.5,
        minimal_dropout: Tunable[float] = 0.01,
        zero_inflation: Tunable[
            Literal["gene", "gene-batch", "gene-label", "gene-cell"]
        ] = "gene",
        **kwargs,
    ) -> None:
        if "reconstruction_loss" in kwargs:
            raise ValueError(
                "No reconstruction loss must be specified for AutoZI : it is 'autozinb'."
            )

        super().__init__(n_input, **kwargs)
        self.zero_inflation = zero_inflation
        self.reconstruction_loss = "autozinb"
        self.minimal_dropout = minimal_dropout

        # Parameters of prior Bernoulli Beta distribution : alpha + beta = 1 if only one is specified
        if beta_prior is None and alpha_prior is not None:
            beta_prior = 1.0 - alpha_prior
        if alpha_prior is None and beta_prior is not None:
            alpha_prior = 1.0 - beta_prior

        # Create parameters for Bernoulli Beta prior and posterior distributions
        # Each parameter, whose values are in (0,1), is encoded as its logit, in the set of real numbers

        if self.zero_inflation == "gene":
            self.alpha_posterior_logit = torch.nn.Parameter(torch.randn(n_input))
            self.beta_posterior_logit = torch.nn.Parameter(torch.randn(n_input))
            if alpha_prior is None:
                self.alpha_prior_logit = torch.nn.Parameter(torch.randn(1))
            else:
                self.register_buffer(
                    "alpha_prior_logit", torch.tensor([logit(alpha_prior)])
                )
            if beta_prior is None:
                self.beta_prior_logit = torch.nn.Parameter(torch.randn(1))
            else:
                self.register_buffer(
                    "beta_prior_logit", torch.tensor([logit(alpha_prior)])
                )

        elif self.zero_inflation == "gene-batch":
            self.alpha_posterior_logit = torch.nn.Parameter(
                torch.randn(n_input, self.n_batch)
            )
            self.beta_posterior_logit = torch.nn.Parameter(
                torch.randn(n_input, self.n_batch)
            )
            if alpha_prior is None:
                self.alpha_prior_logit = torch.nn.parameter(
                    torch.randn(1, self.n_batch)
                )
            else:
                self.register_buffer(
                    "alpha_prior_logit", torch.tensor([logit(alpha_prior)])
                )
            if beta_prior is None:
                self.beta_prior_logit = torch.nn.parameter(torch.randn(1, self.n_batch))
            else:
                self.register_buffer(
                    "beta_prior_logit", torch.tensor([logit(beta_prior)])
                )

        elif self.zero_inflation == "gene-label":
            self.alpha_posterior_logit = torch.nn.Parameter(
                torch.randn(n_input, self.n_labels)
            )
            self.beta_posterior_logit = torch.nn.Parameter(
                torch.randn(n_input, self.n_labels)
            )
            if alpha_prior is None:
                self.alpha_prior_logit = torch.nn.parameter(
                    torch.randn(1, self.n_labels)
                )
            else:
                self.register_buffer(
                    "alpha_prior_logit", torch.tensor([logit(alpha_prior)])
                )
            if beta_prior is None:
                self.beta_prior_logit = torch.nn.parameter(
                    torch.randn(1, self.n_labels)
                )
            else:
                self.register_buffer(
                    "beta_prior_logit", torch.tensor([logit(beta_prior)])
                )

        else:  # gene-cell
            raise Exception("Gene-cell not implemented yet for AutoZI")

    def get_alphas_betas(
        self, as_numpy: bool = True
    ) -> Dict[str, Union[torch.Tensor, np.ndarray]]:
        """Get the parameters of the Bernoulli beta prior and posterior distributions."""
        # Return parameters of Bernoulli Beta distributions in a dictionary
        outputs = {}
        outputs["alpha_posterior"] = torch.sigmoid(self.alpha_posterior_logit)
        outputs["beta_posterior"] = torch.sigmoid(self.beta_posterior_logit)
        outputs["alpha_prior"] = torch.sigmoid(self.alpha_prior_logit)
        outputs["beta_prior"] = torch.sigmoid(self.beta_prior_logit)

        if as_numpy:
            for key, value in outputs.items():
                outputs[key] = (
                    value.detach().cpu().numpy()
                    if value.requires_grad
                    else value.cpu().numpy()
                )

        return outputs

    def sample_from_beta_distribution(
        self,
        alpha: torch.Tensor,
        beta: torch.Tensor,
        eps_gamma: float = 1e-30,
        eps_sample: float = 1e-7,
    ) -> torch.Tensor:
        """Sample from a beta distribution."""
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
        sample_xplusy_log = sample_xy_log_max + torch.log(
            torch.exp(sample_x_log - sample_xy_log_max)
            + torch.exp(sample_y_log - sample_xy_log_max)
        )
        sample_log = sample_x_log - sample_xplusy_log
        sample = eps_sample + (1 - 2 * eps_sample) * torch.exp(sample_log)

        return sample

    def reshape_bernoulli(
        self,
        bernoulli_params: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        y: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Reshape Bernoulli parameters to match the input tensor."""
        if self.zero_inflation == "gene-label":
            one_hot_label = one_hot(y, self.n_labels)
            # If we sampled several random Bernoulli parameters
            if len(bernoulli_params.shape) == 2:
                bernoulli_params = F.linear(one_hot_label, bernoulli_params)
            else:
                bernoulli_params_res = []
                for sample in range(bernoulli_params.shape[0]):
                    bernoulli_params_res.append(
                        F.linear(one_hot_label, bernoulli_params[sample])
                    )
                bernoulli_params = torch.stack(bernoulli_params_res)
        elif self.zero_inflation == "gene-batch":
            one_hot_batch = one_hot(batch_index, self.n_batch)
            if len(bernoulli_params.shape) == 2:
                bernoulli_params = F.linear(one_hot_batch, bernoulli_params)
            # If we sampled several random Bernoulli parameters
            else:
                bernoulli_params_res = []
                for sample in range(bernoulli_params.shape[0]):
                    bernoulli_params_res.append(
                        F.linear(one_hot_batch, bernoulli_params[sample])
                    )
                bernoulli_params = torch.stack(bernoulli_params_res)

        return bernoulli_params

    def sample_bernoulli_params(
        self,
        batch_index: Optional[torch.Tensor] = None,
        y: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> torch.Tensor:
        """Sample Bernoulli parameters from the posterior distribution."""
        outputs = self.get_alphas_betas(as_numpy=False)
        alpha_posterior = outputs["alpha_posterior"]
        beta_posterior = outputs["beta_posterior"]

        if n_samples > 1:
            alpha_posterior = (
                alpha_posterior.unsqueeze(0).expand(
                    (n_samples, alpha_posterior.size(0))
                )
                if self.zero_inflation == "gene"
                else alpha_posterior.unsqueeze(0).expand(
                    (n_samples, alpha_posterior.size(0), alpha_posterior.size(1))
                )
            )
            beta_posterior = (
                beta_posterior.unsqueeze(0).expand((n_samples, beta_posterior.size(0)))
                if self.zero_inflation == "gene"
                else beta_posterior.unsqueeze(0).expand(
                    (n_samples, beta_posterior.size(0), beta_posterior.size(1))
                )
            )

        bernoulli_params = self.sample_from_beta_distribution(
            alpha_posterior, beta_posterior
        )
        bernoulli_params = self.reshape_bernoulli(bernoulli_params, batch_index, y)

        return bernoulli_params

    def rescale_dropout(
        self, px_dropout: torch.Tensor, eps_log: float = 1e-8
    ) -> torch.Tensor:
        """Rescale dropout rate."""
        if self.minimal_dropout > 0.0:
            dropout_prob_rescaled = self.minimal_dropout + (
                1.0 - self.minimal_dropout
            ) * torch.sigmoid(px_dropout)
            px_dropout_rescaled = torch.log(
                dropout_prob_rescaled / (1.0 - dropout_prob_rescaled + eps_log)
            )
        else:
            px_dropout_rescaled = px_dropout
        return px_dropout_rescaled

    def generative(
        self,
        z,
        library,
        batch_index: Optional[torch.Tensor] = None,
        y: Optional[torch.Tensor] = None,
        size_factor=None,
        cont_covs=None,
        cat_covs=None,
        n_samples: int = 1,
        eps_log: float = 1e-8,
    ) -> Dict[str, torch.Tensor]:
        """Run the generative model."""
        outputs = super().generative(
            z=z,
            library=library,
            batch_index=batch_index,
            cont_covs=cont_covs,
            cat_covs=cat_covs,
            y=y,
            size_factor=size_factor,
        )
        # Rescale dropout
        rescaled_dropout = self.rescale_dropout(
            outputs["px"].zi_logits, eps_log=eps_log
        )
        outputs["px"] = ZeroInflatedNegativeBinomial(
            mu=outputs["px"].mu,
            theta=outputs["px"].theta,
            zi_logits=rescaled_dropout,
            scale=outputs["px"].scale,
        )

        # Bernoulli parameters
        outputs["bernoulli_params"] = self.sample_bernoulli_params(
            batch_index, y, n_samples=n_samples
        )
        return outputs

    def compute_global_kl_divergence(self) -> torch.Tensor:
        """Compute global KL divergence."""
        outputs = self.get_alphas_betas(as_numpy=False)
        alpha_posterior = outputs["alpha_posterior"]
        beta_posterior = outputs["beta_posterior"]
        alpha_prior = outputs["alpha_prior"]
        beta_prior = outputs["beta_prior"]

        return kl(
            Beta(alpha_posterior, beta_posterior), Beta(alpha_prior, beta_prior)
        ).sum()

    def get_reconstruction_loss(
        self,
        x: torch.Tensor,
        px_rate: torch.Tensor,
        px_r: torch.Tensor,
        px_dropout: torch.Tensor,
        bernoulli_params: torch.Tensor,
        eps_log: float = 1e-8,
        **kwargs,
    ) -> torch.Tensor:
        """Compute the reconstruction loss."""
        # LLs for NB and ZINB
        ll_zinb = torch.log(
            1.0 - bernoulli_params + eps_log
        ) + ZeroInflatedNegativeBinomial(
            mu=px_rate, theta=px_r, zi_logits=px_dropout
        ).log_prob(
            x
        )
        ll_nb = torch.log(bernoulli_params + eps_log) + NegativeBinomial(
            mu=px_rate, theta=px_r
        ).log_prob(x)

        # Reconstruction loss using a logsumexp-type computation
        ll_max = torch.max(ll_zinb, ll_nb)
        ll_tot = ll_max + torch.log(
            torch.exp(ll_nb - ll_max) + torch.exp(ll_zinb - ll_max)
        )
        reconst_loss = -ll_tot.sum(dim=-1)

        return reconst_loss

    @auto_move_data
    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: int = 1.0,
        n_obs: int = 1.0,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Compute the loss."""
        # Parameters for z latent distribution
        qz = inference_outputs["qz"]
        px = generative_outputs["px"]
        px_rate = px.mu
        px_r = px.theta
        px_dropout = px.zi_logits
        bernoulli_params = generative_outputs["bernoulli_params"]
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        # KL divergences wrt z_n,l_n
        mean = torch.zeros_like(qz.loc)
        scale = torch.ones_like(qz.scale)

        kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)
        if not self.use_observed_lib_size:
            ql = inference_outputs["ql"]

            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)

            kl_divergence_l = kl(
                ql,
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),
            ).sum(dim=1)
        else:
            kl_divergence_l = 0.0

        # KL divergence wrt Bernoulli parameters
        kl_divergence_bernoulli = self.compute_global_kl_divergence()

        # Reconstruction loss
        reconst_loss = self.get_reconstruction_loss(
            x, px_rate, px_r, px_dropout, bernoulli_params
        )

        kl_global = kl_divergence_bernoulli
        kl_local_for_warmup = kl_divergence_z
        kl_local_no_warmup = kl_divergence_l

        weighted_kl_local = kl_weight * kl_local_for_warmup + kl_local_no_warmup
        loss = n_obs * torch.mean(reconst_loss + weighted_kl_local) + kl_global
        kl_local = dict(
            kl_divergence_l=kl_divergence_l, kl_divergence_z=kl_divergence_z
        )
        return LossRecorder(loss, reconst_loss, kl_local, kl_global)
