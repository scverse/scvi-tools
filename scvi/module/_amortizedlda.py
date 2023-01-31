import math
from typing import Dict, Iterable, Optional, Sequence, Tuple, Union

import pyro
import pyro.distributions as dist
import torch
import torch.nn.functional as F
from pyro import poutine
from pyro.infer import Trace_ELBO
from pyro.nn import PyroModule

from scvi._constants import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.module.base import PyroBaseModuleClass, auto_move_data
from scvi.nn import Encoder

_AMORTIZED_LDA_PYRO_MODULE_NAME = "amortized_lda"


class CategoricalBoW(dist.Multinomial):
    """Categorical BoW."""

    def log_prob(self, value):
        """Log probability."""
        if self._validate_args:
            self._validate_sample(value)
        logits, value = dist.util.broadcast_all(self.logits, value)
        logits = logits.clone(memory_format=torch.contiguous_format)
        logits[(value == 0) & (logits == -math.inf)] = 0
        log_powers = (logits * value).sum(-1)
        return log_powers


def logistic_normal_approximation(
    alpha: torch.Tensor,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Returns the mean and standard deviation of the Logistic Normal approximation to the Dirichlet.

    Uses the Laplace approximation of the Logistic Normal distribution to the Dirichlet distribution
    as described in Srivastava et al. https://arxiv.org/pdf/1703.01488.pdf.
    """
    K = alpha.shape[0]
    mu = torch.log(alpha) - torch.log(alpha).sum() / K
    sigma = torch.sqrt((1 - 2 / K) / alpha + torch.sum(1 / alpha) / K**2)
    return mu, sigma


class AmortizedLDAPyroModel(PyroModule):
    """
    A PyroModule that serves as the model for the AmortizedLDAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input features.
    n_topics
        Number of topics/topics to model.
    cell_topic_prior
        Prior of cell topic distribution.
    topic_feature_prior
        Prior of topic feature distribution.
    """

    def __init__(
        self,
        n_input: int,
        n_topics: Tunable[int],
        cell_topic_prior: torch.Tensor,
        topic_feature_prior: torch.Tensor,
    ):
        super().__init__(_AMORTIZED_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_topics = n_topics
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        cell_topic_prior_mu, cell_topic_prior_sigma = logistic_normal_approximation(
            cell_topic_prior
        )
        self.register_buffer(
            "cell_topic_prior_mu",
            cell_topic_prior_mu,
        )
        self.register_buffer(
            "cell_topic_prior_sigma",
            cell_topic_prior_sigma,
        )
        (
            topic_feature_prior_mu,
            topic_feature_prior_sigma,
        ) = logistic_normal_approximation(topic_feature_prior)
        self.register_buffer("topic_feature_prior_mu", topic_feature_prior_mu)
        self.register_buffer("topic_feature_prior_sigma", topic_feature_prior_sigma)

        # Hack: to allow auto_move_data to infer device.
        self._dummy = torch.nn.Parameter(torch.zeros(1), requires_grad=False)

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:

        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        library = torch.sum(x, dim=1)
        return (x, library), {}

    @auto_move_data
    def forward(
        self,
        x: torch.Tensor,
        library: torch.Tensor,
        n_obs: Optional[int] = None,
        kl_weight: float = 1.0,
    ):
        """Forward pass."""
        # Topic feature distributions.
        with pyro.plate("topics", self.n_topics), poutine.scale(None, kl_weight):
            log_topic_feature_dist = pyro.sample(
                "log_topic_feature_dist",
                dist.Normal(
                    self.topic_feature_prior_mu, self.topic_feature_prior_sigma
                ).to_event(1),
            )
            topic_feature_dist = F.softmax(log_topic_feature_dist, dim=1)

        # Cell counts generation.
        max_library_size = int(torch.max(library).item())
        with pyro.plate("cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]):
            with poutine.scale(None, kl_weight):
                log_cell_topic_dist = pyro.sample(
                    "log_cell_topic_dist",
                    dist.Normal(
                        self.cell_topic_prior_mu, self.cell_topic_prior_sigma
                    ).to_event(1),
                )
            cell_topic_dist = F.softmax(log_cell_topic_dist, dim=1)

            pyro.sample(
                "feature_counts",
                CategoricalBoW(max_library_size, cell_topic_dist @ topic_feature_dist),
                obs=x,
            )


class AmortizedLDAPyroGuide(PyroModule):
    """
    A PyroModule that serves as the guide for the AmortizedLDAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input features.
    n_topics
        Number of topics/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    """

    def __init__(self, n_input: int, n_topics: int, n_hidden: int):
        super().__init__(_AMORTIZED_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_topics = n_topics
        self.n_hidden = n_hidden
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        self.encoder = Encoder(n_input, n_topics, distribution="ln", return_dist=True)
        (
            topic_feature_posterior_mu,
            topic_feature_posterior_sigma,
        ) = logistic_normal_approximation(torch.ones(self.n_input))
        self.topic_feature_posterior_mu = torch.nn.Parameter(
            topic_feature_posterior_mu.repeat(self.n_topics, 1)
        )
        self.unconstrained_topic_feature_posterior_sigma = torch.nn.Parameter(
            topic_feature_posterior_sigma.repeat(self.n_topics, 1)
        )

    @property
    def topic_feature_posterior_sigma(self):  # noqa: D102
        return F.softplus(self.unconstrained_topic_feature_posterior_sigma)

    @auto_move_data
    def forward(
        self,
        x: torch.Tensor,
        _library: torch.Tensor,
        n_obs: Optional[int] = None,
        kl_weight: float = 1.0,
    ):
        """Forward pass."""
        # Topic feature distributions.
        with pyro.plate("topics", self.n_topics), poutine.scale(None, kl_weight):
            pyro.sample(
                "log_topic_feature_dist",
                dist.Normal(
                    self.topic_feature_posterior_mu,
                    self.topic_feature_posterior_sigma,
                ).to_event(1),
            )

        # Cell topic distributions guide.
        with pyro.plate(
            "cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]
        ), poutine.scale(None, kl_weight):
            cell_topic_posterior, _ = self.encoder(x)
            cell_topic_posterior_mu = cell_topic_posterior.loc
            cell_topic_posterior_sigma = cell_topic_posterior.scale**2
            pyro.sample(
                "log_cell_topic_dist",
                dist.Normal(
                    cell_topic_posterior_mu, cell_topic_posterior_sigma
                ).to_event(1),
            )


class AmortizedLDAPyroModule(PyroBaseModuleClass):
    """
    An amortized implementation of Latent Dirichlet Allocation :cite:p:`Blei03` implemented in Pyro.

    This module uses auto encoding variational Bayes to optimize the latent variables in the model.
    In particular, a fully-connected neural network is used as an encoder, which takes in feature counts
    as input and outputs the parameters of cell topic distribution. To employ the reparametrization trick
    stably, the Dirichlet priors are approximated by a Logistic-Normal distribution.
    The input feature counts tensor is a cell by features Bag-of-Words(BoW) representation
    of the counts. I.e. the model treats each cell's feature vector as ordered, not
    as unordered as in a Multinomial distribution.

    Parameters
    ----------
    n_input
        Number of input features.
    n_topics
        Number of topics/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_topic_prior
        Prior of cell topic distribution. If `None`, defaults to `1 / n_topics`.
    topic_feature_prior
        Prior of topic feature distribution. If `None`, defaults to `1 / n_topics`.
    """

    def __init__(
        self,
        n_input: int,
        n_topics: int,
        n_hidden: int,
        cell_topic_prior: Optional[Union[float, Sequence[float]]] = None,
        topic_feature_prior: Optional[Union[float, Sequence[float]]] = None,
    ):
        super().__init__()

        self.n_input = n_input
        self.n_topics = n_topics
        self.n_hidden = n_hidden

        if cell_topic_prior is None:
            self.cell_topic_prior = torch.full((n_topics,), 1 / self.n_topics)
        elif isinstance(cell_topic_prior, float):
            self.cell_topic_prior = torch.full((n_topics,), cell_topic_prior)
        else:
            self.cell_topic_prior = torch.tensor(cell_topic_prior)

        if topic_feature_prior is None:
            self.topic_feature_prior = torch.full((n_input,), 1 / self.n_topics)
        elif isinstance(topic_feature_prior, float):
            self.topic_feature_prior = torch.full((n_input,), topic_feature_prior)
        else:
            self.topic_feature_prior = torch.tensor(topic_feature_prior)

        self._model = AmortizedLDAPyroModel(
            self.n_input,
            self.n_topics,
            self.cell_topic_prior,
            self.topic_feature_prior,
        )
        self._guide = AmortizedLDAPyroGuide(self.n_input, self.n_topics, self.n_hidden)
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):  # noqa: D102
        return self._model

    @property
    def guide(self):  # noqa: D102
        return self._guide

    def topic_by_feature(self, n_samples: int) -> torch.Tensor:
        """
        Gets a Monte-Carlo estimate of the expectation of the topic by feature matrix.

        Assumes the module has already been trained.

        Parameters
        ----------
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.

        Returns
        -------
        A `n_topics x n_input` tensor containing the topic by feature matrix.
        """
        topic_feature_posterior_mu, topic_feature_posterior_sigma = (
            self.guide.topic_feature_posterior_mu.detach().cpu(),
            self.guide.topic_feature_posterior_sigma.detach().cpu(),
        )
        return torch.mean(
            F.softmax(
                dist.Normal(
                    topic_feature_posterior_mu,
                    topic_feature_posterior_sigma,
                ).sample(sample_shape=torch.Size((n_samples,))),
                dim=2,
            ),
            dim=0,
        )

    @auto_move_data
    @torch.inference_mode()
    def get_topic_distribution(self, x: torch.Tensor, n_samples: int) -> torch.Tensor:
        """
        Converts `x` to its inferred topic distribution.

        Parameters
        ----------
        x
            Counts tensor.
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.

        Returns
        -------
        A `x.shape[0] x n_topics` tensor containing the normalized topic distribution.
        """
        cell_topic_dist, _ = self.guide.encoder(x)
        cell_topic_dist_mu = cell_topic_dist.loc.detach().cpu()
        cell_topic_dist_sigma = 2.0 * cell_topic_dist.scale.log()
        cell_topic_dist_sigma = F.softplus(cell_topic_dist_sigma.detach().cpu())
        return torch.mean(
            F.softmax(
                dist.Normal(cell_topic_dist_mu, cell_topic_dist_sigma).sample(
                    sample_shape=torch.Size((n_samples,))
                ),
                dim=2,
            ),
            dim=0,
        )

    @auto_move_data
    @torch.inference_mode()
    def get_elbo(self, x: torch.Tensor, library: torch.Tensor, n_obs: int) -> float:
        """
        Computes ELBO.

        Parameters
        ----------
        x
            Counts tensor.
        library
            Library sizes for each cell.
        n_obs
            Size of full batch. If n_obs < x.shape[0], ELBO is scaled by (n_obs / x.shape[0]).

        Returns
        -------
        The positive ELBO.
        """
        return Trace_ELBO().loss(self.model, self.guide, x, library, n_obs=n_obs)
