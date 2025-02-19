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

    def kl(self, qz: torch.Tensor, z: None = None) -> torch.Tensor:
        """Compute KL div between std. normal prior and the posterior distn.

        Parameters
        ----------
        qz
            Posterior distribution.
        z
            Ignored.

        Returns
        -------
        KL divergence.
        """
        # 1 x N
        return kl_divergence(qz, Normal(torch.zeros_like(qz.loc), torch.ones_like(qz.loc))).sum(
            dim=1
        )


class VampPrior(Prior):
    """VampPrior.

    Adapted from  a
    `blog post
    <https://github.com/jmtomczak/intro_dgm/blob/main/vaes/vae_priors_example.ipynb>`_
    of the original VampPrior author.

    Parameters
    ----------
    n_components
        Number of prior components.
    encoder
        The encoder.
    x
        Expression data for pseudoinputs initialisation.
    n_cat_list
        The number of categorical covariates and
        the number of category levels.
        A list containing, for each covariate of interest,
        the number of categories.
    batch_index
        Batch index for pseudoinputs initialisation.
    cat
        List of categorical covariates for pseudoinputs initialisation.
        Includes all covariates that will be one-hot encoded by the ``encoder``,
        including the batch.
    cont
        Continuous covariates for pseudoinputs initialisation.
    trainable_priors
        Are pseudoinput parameters trainable or fixed.
    """

    # K - components, I - inputs, L - latent, N - samples

    def __init__(
        self,
        n_components: int,
        encoder: torch.nn.Module,
        x: torch.Tensor,
        n_cat_list: list[int],
        batch_index: torch.Tensor,
        cat_list: list[torch.Tensor],
        cont: torch.Tensor | None = None,
        trainable_priors: bool = True,
    ):
        super().__init__()

        self.cat_list = []
        self.encoder = encoder

        # Make pseudoinputs into parameters
        # X
        assert n_components == x.shape[0]
        self.u = torch.nn.Parameter(x, requires_grad=trainable_priors)  # K x I
        # Cat
        cat_list = [batch_index] + cat_list
        assert all(cat.shape[0] == n_components for cat in cat_list)
        # For categorical covariates, since scvi-tools one-hot encodes
        # them in the layers, we need to create a multinomial distribution
        # from which we can sample categories for layers input
        # Initialise the multinomial distribution weights based on
        # one-hot encoding of pseudoinput categories
        self.u_cat = torch.nn.ParameterList(
            [
                torch.nn.Parameter(
                    torch.nn.functional.one_hot(cat.squeeze(-1), n).float(),
                    # K x C_cat_onehot
                    requires_grad=trainable_priors,
                )
                for cat, n in zip(cat_list, n_cat_list, strict=False)
                # K x C_cat
            ]
        )
        # Cont
        if cont is None:
            self.u_cont = None
        else:
            assert n_components == cont.shape[0]
            self.u_cont = torch.nn.Parameter(cont, requires_grad=trainable_priors)  # K x C_cont

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
        batch_index, cat_list = cat_list[0], cat_list[1:]
        z = self.encoder(x=self.u, batch_index=batch_index, cat_list=cat_list, cont=self.u_cont)
        self.encoder.train(original_mode)
        return z["q_dist"].loc, z["q_dist"].scale  # (K x L), (K x L)

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

    def kl(self, qz: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        """Compute KL div. between VampPrior and the posterior distribution.

        Parameters
        ----------
        qz
            Posterior distribution.
        z
            Sample from the posterior distribution.

        Returns
        -------
        KL divergence.
        """
        return (qz.log_prob(z) - self.log_prob(z)).sum(1)
