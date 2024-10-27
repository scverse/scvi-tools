from __future__ import annotations

import abc
from abc import abstractmethod

import torch
from torch.distributions import Normal, kl_divergence


class Prior(torch.nn.Module, abc.ABC):
    """Abstract class for prior distributions."""

    @abstractmethod
    def kl(
        self,
        m_q: torch.Tensor,
        v_q: torch.Tensor,
        z: torch.Tensor,
    ) -> torch.Tensor:
        """Compute KL divergence between prior and posterior distribution.

        Parameters
        ----------
        m_q
            Posterior distribution mean.
        v_q
            Posterior distribution variance.
        z
            Sample from the posterior distribution.

        Returns
        -------
        KL divergence.
        """
        pass


class StandardPrior(Prior):
    """Standard prior distribution."""

    def kl(self, m_q: torch.Tensor, v_q: torch.Tensor, z: None = None) -> torch.Tensor:
        """Compute KL divergence between standard normal prior and the posterior distribution.

        Parameters
        ----------
        m_q
            Posterior distribution mean.
        v_q
            Posterior distribution variance.
        z
            Ignored.

        Returns
        -------
        KL divergence.
        """
        # 1 x N
        return kl_divergence(
            Normal(m_q, v_q.sqrt()), Normal(torch.zeros_like(m_q), torch.ones_like(v_q))
        ).sum(dim=1)


class VampPrior(Prior):
    """VampPrior.

    Adapted from  a
    `blog post <https://github.com/jmtomczak/intro_dgm/blob/main/vaes/vae_priors_example.ipynb>`_
    of the original VampPrior author.

    Parameters
    ----------
    n_components
        Number of prior components.
    encoder
        The encoder.
    data_x
        Expression data for pseudoinputs initialisation.
    n_cat_list
        The number of categorical covariates and
        the number of category levels.
        A list containing, for each covariate of interest,
        the number of categories.
    data_cat
        List of categorical covariates for pseudoinputs initialisation.
        Includes all covariates that will be one-hot encoded by the ``encoder``,
        including the batch.
    data_cont
        Continuous covariates for pseudoinputs initialisation.
    trainable_priors
        Are pseudoinput parameters trainable or fixed.
    """

    # K - components, I - inputs, L - latent, N - samples

    def __init__(
        self,
        n_components: int,
        encoder: torch.nn.Module,
        data_x: torch.Tensor,
        n_cat_list: list[int],
        data_cat: list[torch.Tensor],
        data_cont: torch.Tensor | None = None,
        trainable_priors: bool = True,
    ):
        super().__init__()

        self.encoder = encoder

        # Make pseudoinputs into parameters
        # X
        assert n_components == data_x.shape[0]
        self.u = torch.nn.Parameter(data_x, requires_grad=trainable_priors)  # K x I
        # Cat
        assert all([cat.shape[0] == n_components for cat in data_cat])
        # For categorical covariates, since scvi-tools one-hot encodes
        # them in the layers, we need to create a multinomial distn
        # from which we can sample categories for layers input
        # Initialise the multinomial distn weights based on
        # one-hot encoding of pseudoinput categories
        self.u_cat = torch.nn.ParameterList(
            [
                torch.nn.Parameter(
                    torch.nn.functional.one_hot(cat.squeeze(-1), n).float(),  # K x C_cat_onehot
                    requires_grad=trainable_priors,
                )
                for cat, n in zip(data_cat, n_cat_list, strict=False)  # K x C_cat
            ]
        )
        # Cont
        if data_cont is None:
            self.u_cont = None
        else:
            assert n_components == data_cont.shape[0]
            self.u_cont = torch.nn.Parameter(
                data_cont, requires_grad=trainable_priors
            )  # K x C_cont

        # mixing weights
        self.w = torch.nn.Parameter(torch.zeros(n_components, 1, 1))  # K x 1 x 1

    def get_params(self) -> tuple[torch.Tensor, torch.Tensor]:
        """Get posterior of pseudoinputs.

        Returns
        -------
        Posterior representation mean and variance for each pseudoinput.
        """
        # u, u_cov -> encoder -> mean, var
        original_mode = self.encoder.training
        self.encoder.train(False)
        # Convert category weights to categories
        cat_list = [torch.multinomial(cat, num_samples=1) for cat in self.u_cat]
        z = self.encoder(x=self.u, cat_list=cat_list, cont=self.u_cont)
        self.encoder.train(original_mode)
        return z["y_m"], z["y_v"]  # (K x L), (K x L)

    def log_prob(self, z: torch.Tensor) -> torch.Tensor:
        """Log probability of posterior sample under the prior.

        Parameters
        ----------
        z
            Latent embedding of samples.

        Returns
        -------
            Log probability of every sample.
            dim = n_samples * n_latent_dimensions
        """
        # Mixture of gaussian computed on K x N x L
        z = z.unsqueeze(0)  # 1 x N x L

        # Get pseudoinputs posteriors which are prior params
        m_p, v_p = self.get_params()  # (K x L), (K x L)
        m_p = m_p.unsqueeze(1)  # K x 1 x L
        v_p = v_p.unsqueeze(1)  # K x 1 x L

        # mixing probabilities
        w = torch.nn.functional.softmax(self.w, dim=0)  # K x 1 x 1

        # sum of log_p across components weighted by w
        log_prob = Normal(m_p, v_p.sqrt()).log_prob(z) + torch.log(w)  # K x N x L
        log_prob = torch.logsumexp(log_prob, dim=0, keepdim=False)  # N x L

        return log_prob  # N x L

    def kl(self, m_q: torch.Tensor, v_q: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        """Compute KL divergence between VampPrior and the posterior distribution.

        Parameters
        ----------
        m_q
            Posterior distribution mean.
        v_q
            Posterior distribution variance.
        z
            Sample from the posterior distribution.

        Returns
        -------
        KL divergence.
        """
        return (Normal(m_q, v_q.sqrt()).log_prob(z) - self.log_prob(z)).sum(1)
