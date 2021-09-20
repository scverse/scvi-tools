import math
from typing import Dict, Iterable, Optional, Sequence, Tuple, Union

import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
import torch
import torch.nn.functional as F
from pyro.infer import Trace_ELBO
from pyro.nn import PyroModule

from scvi._constants import _CONSTANTS
from scvi.module.base import PyroBaseModuleClass, auto_move_data
from scvi.nn import Encoder

_LDA_PYRO_MODULE_NAME = "lda"


class CategoricalBoW(dist.Multinomial):
    def log_prob(self, value):
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
    K = alpha.shape[0]
    mu = torch.log(alpha) - torch.log(alpha).sum() / K
    sigma = (1 - 2 / K) / alpha + torch.sum(1 / alpha) / K ** 2
    return mu, sigma


class LDAPyroModel(PyroModule):
    """
    A PyroModule that serves as the model for the LDAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_topics
        Number of topics/topics to model.
    cell_topic_prior
        Prior of cell topic distribution.
    topic_gene_prior
        Prior of topic gene distribution.
    """

    def __init__(
        self,
        n_input: int,
        n_topics: int,
        cell_topic_prior: torch.Tensor,
        topic_gene_prior: torch.Tensor,
    ):
        super().__init__(_LDA_PYRO_MODULE_NAME)

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
        topic_gene_prior_mu, topic_gene_prior_sigma = logistic_normal_approximation(
            topic_gene_prior
        )
        self.register_buffer("topic_gene_prior_mu", topic_gene_prior_mu)
        self.register_buffer("topic_gene_prior_sigma", topic_gene_prior_sigma)

        # Hack: to allow auto_move_data to infer device.
        self._dummy = torch.nn.Parameter(torch.zeros(1), requires_grad=False)

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:

        x = tensor_dict[_CONSTANTS.X_KEY]
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
        # Topic gene distributions.
        with pyro.plate("topics", self.n_topics), poutine.scale(None, kl_weight):
            log_topic_gene_dist = pyro.sample(
                "log_topic_gene_dist",
                dist.Normal(
                    self.topic_gene_prior_mu, self.topic_gene_prior_sigma
                ).to_event(1),
            )
            topic_gene_dist = F.softmax(log_topic_gene_dist, dim=1)

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
                "gene_counts",
                CategoricalBoW(max_library_size, cell_topic_dist @ topic_gene_dist),
                obs=x,
            )


class LDAPyroGuide(PyroModule):
    """
    A PyroModule that serves as the guide for the LDAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_topics
        Number of topics/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    """

    def __init__(self, n_input: int, n_topics: int, n_hidden: int):
        super().__init__(_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_topics = n_topics
        self.n_hidden = n_hidden
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        self.encoder = Encoder(n_input, n_topics, distribution="ln")
        (
            topic_gene_posterior_mu,
            topic_gene_posterior_sigma,
        ) = logistic_normal_approximation(torch.ones(self.n_input))
        self.topic_gene_posterior_mu = torch.nn.Parameter(
            topic_gene_posterior_mu.repeat(self.n_topics, 1)
        )
        self.unconstrained_topic_gene_posterior_sigma = torch.nn.Parameter(
            topic_gene_posterior_sigma.repeat(self.n_topics, 1)
        )

    @property
    def topic_gene_posterior_sigma(self):
        return F.softplus(self.unconstrained_topic_gene_posterior_sigma)

    @auto_move_data
    def forward(
        self,
        x: torch.Tensor,
        _library: torch.Tensor,
        n_obs: Optional[int] = None,
        kl_weight: float = 1.0,
    ):
        # Topic gene distributions.
        with pyro.plate("topics", self.n_topics), poutine.scale(None, kl_weight):
            pyro.sample(
                "log_topic_gene_dist",
                dist.Normal(
                    self.topic_gene_posterior_mu,
                    self.topic_gene_posterior_sigma,
                ).to_event(1),
            )

        # Cell topic distributions guide.
        with pyro.plate(
            "cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]
        ), poutine.scale(None, kl_weight):
            cell_topic_posterior_mu, cell_topic_posterior_sigma, _ = self.encoder(x)
            pyro.sample(
                "log_cell_topic_dist",
                dist.Normal(
                    cell_topic_posterior_mu, F.softplus(cell_topic_posterior_sigma)
                ).to_event(1),
            )


class LDAPyroModule(PyroBaseModuleClass):
    """
    Latent Dirichlet Allocation [Blei03]_ implemented in Pyro.

    This model uses auto encoding variational Bayes to optimize the latent variables in the model.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_topics
        Number of topics/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_topic_prior
        Prior of cell topic distribution. If `None`, defaults to `1 / n_topics`.
    topic_gene_prior
        Prior of topic gene distribution. If `None`, defaults to `1 / n_topics`.
    """

    def __init__(
        self,
        n_input: int,
        n_topics: int,
        n_hidden: int,
        cell_topic_prior: Optional[Union[float, Sequence[float]]] = None,
        topic_gene_prior: Optional[Union[float, Sequence[float]]] = None,
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

        if topic_gene_prior is None:
            self.topic_gene_prior = torch.full((n_input,), 1 / self.n_topics)
        elif isinstance(topic_gene_prior, float):
            self.topic_gene_prior = torch.full((n_input,), topic_gene_prior)
        else:
            self.topic_gene_prior = torch.tensor(topic_gene_prior)

        self._model = LDAPyroModel(
            self.n_input,
            self.n_topics,
            self.cell_topic_prior,
            self.topic_gene_prior,
        )
        self._guide = LDAPyroGuide(self.n_input, self.n_topics, self.n_hidden)
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide

    def topic_by_gene(self, give_mean: bool = True) -> torch.Tensor:
        """
        Gets the topic by gene matrix.

        Assumes the module has already been trained.

        Parameters
        ----------
        give_mean
            Give mean of distribution if True or sample from it.

        Returns
        -------
        A `n_topics x n_input` tensor containing the topic by gene matrix.
        """
        topic_gene_posterior_mu, topic_gene_posterior_sigma = (
            self.guide.topic_gene_posterior_mu.detach().cpu(),
            self.guide.topic_gene_posterior_sigma.detach().cpu(),
        )
        if give_mean:
            return F.softmax(topic_gene_posterior_mu, dim=1)
        return F.softmax(
            torch.normal(topic_gene_posterior_mu, topic_gene_posterior_sigma),
            dim=1,
        )

    @auto_move_data
    @torch.no_grad()
    def get_topic_distribution(
        self, x: torch.Tensor, give_mean: bool = True
    ) -> torch.Tensor:
        """
        Converts `x` to its inferred topic distribution.

        Parameters
        ----------
        give_mean
            Give mean of distribution or sample from it.

        Returns
        -------
        A `x.shape[0] x n_topics` tensor containing the normalized topic distribution.
        """
        (
            cell_topic_dist_mu,
            _,
            cell_topic_dist_sample,
        ) = self.guide.encoder(x)
        if give_mean:
            return F.softmax(cell_topic_dist_mu.detach().cpu(), dim=1)
        return cell_topic_dist_sample.detach().cpu()

    @auto_move_data
    @torch.no_grad()
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
