from __future__ import annotations

import abc
from abc import abstractmethod
from typing import TYPE_CHECKING

import torch
from torch.distributions import (
    Categorical,
    Independent,
    MixtureSameFamily,
    Normal,
    kl_divergence,
)

from scvi.module._constants import MODULE_KEYS

if TYPE_CHECKING:
    from collections.abc import Callable

    from torch.distributions import Distribution


class Prior(torch.nn.Module, abc.ABC):
    """Abstract class for prior distributions."""

    @abstractmethod
    def kl(
        self,
        qz: Distribution,
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


class GaussianPrior(Prior):
    """Standard prior distribution."""

    def get_prior(self, z, **kwargs) -> Distribution:
        """Standard normal prior distribution.

        Parameters
        ----------
        z
            Latent embedding of samples.
        """
        return Normal(torch.zeros_like(z), torch.ones_like(z))

    def kl(self, qz: Distribution, **kwargs) -> torch.Tensor:
        """Compute KL div between std. normal prior and the posterior distn.

        Parameters
        ----------
        qz
            Posterior distribution.
        kwargs
            Ignored. Compatibility with other priors.
        """
        prior = self.get_prior(qz.sample())
        return kl_divergence(qz, prior).sum(dim=1)


class MogPrior(Prior):
    """Mixture-of-Gaussian prior.

    Parameters
    ----------
    n_components
        Number of prior components.
    n_latent
        Number of latent dimensions.
    celltype_bias
        Biased logits to match one cell-type to each component.
    """

    def __init__(
        self,
        n_components: int,
        n_latent: int,
        celltype_bias: bool = False,
    ):
        super().__init__()
        self.prior_means = torch.nn.Parameter(0.1 * torch.randn([n_components, n_latent]))
        self.prior_log_scales = torch.nn.Parameter(torch.zeros([n_components, n_latent]) - 1.0)
        self.prior_logits = torch.nn.Parameter(torch.zeros([n_components]))
        self.celltype_bias = celltype_bias
        if celltype_bias:
            self.n_labels = n_components

    def get_prior(self, z: torch.Tensor, y: torch.Tensor | None = None, **kwargs) -> torch.Tensor:
        """Learned prior distribution.

        Parameters
        ----------
        z
            Latent embedding of samples.
        y
            Cell-type labels. By default None. Required if using cell-type bias.
        """
        if self.celltype_bias:
            logits_input = torch.where(
                y < self.n_labels,
                torch.nn.functional.one_hot(y.ravel(), num_classes=self.n_labels),
                torch.zeros(y.shape[0], self.n_labels, device=y.device),
            )
            prior_logits = self.prior_logits + 10 * logits_input
            prior_means = self.prior_means.expand(y.shape[0], -1, -1)
        else:
            prior_logits = self.prior_logits
            prior_means = self.prior_means
        cats = Categorical(logits=prior_logits)
        normal_dists = Independent(
            Normal(prior_means, torch.exp(self.prior_log_scales) + 1e-4),
            reinterpreted_batch_ndims=1,
        )
        prior = MixtureSameFamily(cats, normal_dists)

        return prior

    def kl(self, qz: torch.Tensor, labels: torch.Tensor | None = None, **kwargs) -> torch.Tensor:
        """Compute KL div. between Mixture-of-Gaussian prior and the posterior distribution.

        Parameters
        ----------
        qz
            Posterior distribution.
        labels
            Cell-type labels. By default None. Required if using cell-type bias.
        """
        z = qz.rsample(sample_shape=(30,))
        prior = self.get_prior(z, y=labels)
        return (qz.log_prob(z).sum(-1) - prior.log_prob(z)).mean(0)


class VampPrior(Prior):
    """VampPrior.

    Adapted from a
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
        n_cat_list: list[int],
        inference: Callable,
        encoder: torch.nn.Module,
        pseudoinputs: dict,
        trainable_priors: bool = True,
        additional_categorical_covariates: list[str] | None = None,
    ):
        super().__init__()

        if pseudoinputs["cat_covs"] is None:
            pseudoinputs["cat_covs"] = []

        self.inference = inference
        self.encoder = encoder

        # Make pseudoinputs into parameters
        assert n_components == pseudoinputs["x"].shape[0]
        self.pseudoinputs = pseudoinputs
        self.pseudoinput_x = torch.nn.Parameter(
            torch.log(pseudoinputs.pop("x") + 1e-6), requires_grad=trainable_priors
        )  # K x I
        cat_list = [pseudoinputs["batch_index"]] + pseudoinputs["cat_covs"]
        if additional_categorical_covariates is not None:
            for i in additional_categorical_covariates:
                cat_list.append(pseudoinputs[i])
        else:
            additional_categorical_covariates = []
        self.additional_categorical_covariates = additional_categorical_covariates
        assert all(cat.shape[0] == n_components for cat in cat_list)
        # For categorical covariates, since scvi-tools one-hot encodes
        # them in the layers, we need to create a multinomial distribution
        # from which we can sample categories for layers input
        # Initialise the multinomial distribution weights based on
        # one-hot encoding of pseudoinput categories
        self.u_cat = torch.nn.ParameterList(
            [
                torch.nn.Parameter(
                    torch.nn.functional.one_hot(cat.long().squeeze(-1), n).float(),
                    # K x C_cat_onehot
                    requires_grad=trainable_priors,
                )
                for cat, n in zip(cat_list, n_cat_list, strict=True)
                # K x C_cat
            ]
        )
        # Continuous covariates
        if self.pseudoinputs["cont_covs"] is None:
            self.u_cont = None
        else:
            assert n_components == self.pseudoinputs["cont_covs"].shape[0]
            self.pseudoinputs["cont_covs"] = torch.nn.Parameter(
                self.pseudoinputs["cont_covs"], requires_grad=trainable_priors
            )  # K x C_cont

        # mixing weights
        self.w = torch.nn.Parameter(torch.zeros(n_components))  # K x 1

    def get_prior(self, **kwargs) -> tuple[torch.Tensor, torch.Tensor]:
        """Get prior latent distribution."""
        # u, u_cov -> encoder -> mean, var
        # Convert category weights to categories
        cat_list = [torch.multinomial(cat, num_samples=1) for cat in self.u_cat]
        for key, value in self.pseudoinputs.items():
            if isinstance(value, torch.Tensor):
                if value.device != self.w.device:
                    self.pseudoinputs[key] = value.to(self.w.device)
        for i in self.additional_categorical_covariates:
            self.pseudoinputs[i] = cat_list[-1]
            cat_list = cat_list[:-1]

        if len(cat_list) > 1:
            self.pseudoinputs_["batch_index"], self.pseudoinputs["cat_covs"] = (
                cat_list[0],
                cat_list[1:],
            )
        else:
            self.pseudoinputs["batch_index"], self.pseudoinputs["cat_covs"] = cat_list[0], None

        # Encoder in evaluation mode.
        original_mode = self.encoder.training  # Store original mode
        self.encoder.eval()  # Set to evaluation mode
        try:
            inference_out = self.inference(torch.exp(self.pseudoinput_x), **self.pseudoinputs)
        finally:
            self.encoder.train(original_mode)  # Restore original mode

        cat = Categorical(logits=self.w)  # K x 1
        normal_dists = Independent(
            inference_out[MODULE_KEYS.QZ_KEY],  # K x L
            reinterpreted_batch_ndims=1,
        )
        prior = MixtureSameFamily(cat, normal_dists)

        return prior

    def kl(self, qz: torch.Tensor, z: torch.Tensor, **kwargs) -> torch.Tensor:
        """Compute KL div. between VampPrior and the posterior distribution.

        Parameters
        ----------
        qz
            Posterior distribution.
        z
            Sample from the posterior distribution.
        kwargs
            Ignored, compatibility with other priors.

        Returns
        -------
        KL divergence.
        """
        prior = self.get_prior()
        z = qz.rsample(sample_shape=(30,))
        return (qz.log_prob(z).sum(-1) - prior.log_prob(z)).mean(0)
