import math
from typing import Dict, Iterable, Optional, Sequence, Union

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
import torch.nn.functional as F
from pyro.infer import Trace_ELBO
from pyro.nn import PyroModule

from scvi._constants import _CONSTANTS
from scvi.module.base import PyroBaseModuleClass, auto_move_data
from scvi.nn import FCLayers

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
        cell_topic_prior: Sequence[float],
        topic_gene_prior: Sequence[float],
    ):
        super().__init__(_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_topics = n_topics
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        self.register_buffer(
            "cell_topic_prior",
            torch.FloatTensor(cell_topic_prior),
        )
        self.register_buffer("topic_gene_prior", torch.FloatTensor(topic_gene_prior))

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
        self, x: torch.Tensor, library: torch.Tensor, n_obs: Optional[int] = None
    ):
        # Topic gene distributions.
        with pyro.plate("topics", self.n_topics):
            topic_gene_dist = pyro.sample(
                "topic_gene_dist",
                dist.Dirichlet(torch.clamp(self.topic_gene_prior, min=1e-8)),
            )

        # Cell counts generation.
        max_library_size = int(torch.max(library).item())
        with pyro.plate("cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]):
            cell_topic_dist = pyro.sample(
                "cell_topic_dist",
                dist.Dirichlet(torch.clamp(self.cell_topic_prior, min=1e-8)),
            )

            pyro.sample(
                "gene_counts",
                CategoricalBoW(max_library_size, cell_topic_dist @ topic_gene_dist),
                obs=x,
            )


class CellTopicDistPriorEncoder(nn.Module):
    """
    The neural network encoder used in LDAPyroGuide which outputs cell topic posterior estimate.

    Composed of a single hidden layer fully connected neural network followed by a
    log transformation and softplus.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_topics
        Number of topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    """

    def __init__(self, n_input: int, n_topics: int, n_hidden: int):
        super().__init__()

        self.encoder = FCLayers(
            n_in=n_input,
            n_out=n_topics,
            n_hidden=n_hidden,
            n_layers=2,
            inject_covariates=False,
        )

    @auto_move_data
    def forward(self, x: torch.Tensor):
        return F.softplus(self.encoder(torch.log(1 + x)))


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

        self.encoder = CellTopicDistPriorEncoder(n_input, n_topics, n_hidden)
        self.unconstrained_topic_gene_posterior = torch.nn.Parameter(
            torch.ones(self.n_topics, self.n_input),
        )

    @property
    def topic_gene_posterior(self):
        return F.softplus(self.unconstrained_topic_gene_posterior)

    @auto_move_data
    def forward(
        self, x: torch.Tensor, _library: torch.Tensor, n_obs: Optional[int] = None
    ):
        # Topic gene distributions.
        with pyro.plate("topics", self.n_topics):
            pyro.sample(
                "topic_gene_dist",
                dist.Dirichlet(torch.clamp(self.topic_gene_posterior, min=1e-8)),
            )

        # Cell topic distributions guide.
        with pyro.plate("cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]):
            cell_topic_posterior = self.encoder(x)
            pyro.sample(
                "cell_topic_dist",
                dist.Dirichlet(torch.clamp(cell_topic_posterior, min=1e-8)),
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

    @property
    def topic_by_gene(self) -> torch.Tensor:
        """
        Gets the topic by gene matrix.

        Assumes the module has already been trained.

        Returns
        -------
        A `n_topics x n_input` tensor containing the topic by gene matrix.
        """
        return self.guide.topic_gene_posterior.detach().cpu()

    @auto_move_data
    @torch.no_grad()
    def _get_unnormalized_topic_distribution(self, x: torch.Tensor) -> torch.Tensor:
        """
        Converts `x` to its inferred unnormalized topic distribution.

        Returns
        -------
        A `x.shape[0] x n_topics` tensor containing the unnormalized topic distribution.
        """
        return self.guide.encoder(x).detach().cpu()

    def get_topic_distribution(self, x: torch.Tensor) -> torch.Tensor:
        """
        Converts `x` to its inferred normalized topic distribution.

        Returns
        -------
        A `x.shape[0] x n_topics` tensor containing the normalized topic distribution.
        """
        cell_topic_unnormalized_dist = self._get_unnormalized_topic_distribution(x)
        return (
            cell_topic_unnormalized_dist
            / cell_topic_unnormalized_dist.sum(axis=1)[:, np.newaxis]
        )

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
