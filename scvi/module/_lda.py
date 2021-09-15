import math
from typing import Dict, Iterable, Optional, Union

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
    n_components
        Number of components/topics to model.
    cell_component_prior
        Prior of cell component distribution.
    component_gene_prior
        Prior of component gene distribution.
    """

    def __init__(
        self,
        n_input: int,
        n_components: int,
        cell_component_prior: float,
        component_gene_prior: float,
    ):
        super().__init__(_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_components = n_components
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        self.register_buffer(
            "cell_component_prior",
            torch.full((self.n_components,), cell_component_prior),
        )
        self.register_buffer(
            "component_gene_prior", torch.full((self.n_input,), component_gene_prior)
        )

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
        # Component gene distributions.
        with pyro.plate("components", self.n_components):
            component_gene_dist = pyro.sample(
                "component_gene_dist",
                dist.Dirichlet(self.component_gene_prior),
            ).to(x.device)

        # Cell counts generation.
        max_library_size = int(torch.max(library).item())
        with pyro.plate("cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]):
            cell_component_dist = pyro.sample(
                "cell_component_dist", dist.Dirichlet(self.cell_component_prior)
            )

            pyro.sample(
                "gene_counts",
                CategoricalBoW(
                    max_library_size, cell_component_dist @ component_gene_dist
                ),
                obs=x,
            )


class CellComponentDistPriorEncoder(nn.Module):
    """
    The neural network encoder used in LDAPyroGuide which outputs cell component posterior estimate.

    Composed of a single hidden layer fully connected neural network followed by a
    log transformation and softplus.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_components
        Number of components/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    """

    def __init__(self, n_input: int, n_components: int, n_hidden: int):
        super().__init__()

        self.encoder = FCLayers(
            n_in=n_input,
            n_out=n_components,
            n_hidden=n_hidden,
            n_layers=1,
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
    n_components
        Number of components/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    """

    def __init__(self, n_input: int, n_components: int, n_hidden: int):
        super().__init__(_LDA_PYRO_MODULE_NAME)

        self.n_input = n_input
        self.n_components = n_components
        self.n_hidden = n_hidden
        # Populated by PyroTrainingPlan.
        self.n_obs = None

        self.encoder = CellComponentDistPriorEncoder(n_input, n_components, n_hidden)
        self.unconstrained_component_gene_posterior = torch.nn.Parameter(
            torch.ones(self.n_components, self.n_input),
        )

    @property
    def component_gene_posterior(self):
        return F.softplus(self.unconstrained_component_gene_posterior)

    @auto_move_data
    def forward(
        self, x: torch.Tensor, _library: torch.Tensor, n_obs: Optional[int] = None
    ):
        # Component gene distributions.
        with pyro.plate("components", self.n_components):
            pyro.sample(
                "component_gene_dist",
                dist.Dirichlet(self.component_gene_posterior),
            )

        # Cell component distributions guide.
        with pyro.plate("cells", size=n_obs or self.n_obs, subsample_size=x.shape[0]):
            cell_component_posterior = self.encoder(x)
            pyro.sample("cell_component_dist", dist.Dirichlet(cell_component_posterior))


class LDAPyroModule(PyroBaseModuleClass):
    """
    Latent Dirichlet Allocation [Blei03]_ implemented in Pyro.

    This model uses auto encoding variational Bayes to optimize the latent variables in the model.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_components
        Number of components/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_component_prior
        Prior of cell component distribution. If `None`, defaults to `1 / n_components`.
    component_gene_prior
        Prior of component gene distribution. If `None`, defaults to `1 / n_components`.
    """

    def __init__(
        self,
        n_input: int,
        n_components: int,
        n_hidden: int,
        cell_component_prior: Optional[float] = None,
        component_gene_prior: Optional[float] = None,
    ):
        super().__init__()

        self.n_input = n_input
        self.n_components = n_components
        self.n_hidden = n_hidden

        self.cell_component_prior = cell_component_prior or 1 / self.n_components
        self.component_gene_prior = component_gene_prior or 1 / self.n_components
        self._model = LDAPyroModel(
            self.n_input,
            self.n_components,
            self.cell_component_prior,
            self.component_gene_prior,
        )
        self._guide = LDAPyroGuide(self.n_input, self.n_components, self.n_hidden)
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide

    @property
    def components(self) -> torch.Tensor:
        """
        Gets the component to gene transition matrix from the Pyro parameter store.

        Assumes the module has already been trained.

        Returns
        -------
        A `n_components x n_input` tensor containing the component to gene transition matrix.
        """
        return self.guide.component_gene_posterior.detach().cpu()

    @auto_move_data
    @torch.no_grad()
    def _unnormalized_transform(self, x: torch.Tensor) -> torch.Tensor:
        """
        Converts `x` to its inferred unnormalized component distribution.

        Returns
        -------
        A `x.shape[0] x n_components` tensor containing the unnormalized component distribution.
        """
        return self.guide.encoder(x).detach().cpu()

    def transform(self, x: torch.Tensor) -> torch.Tensor:
        """
        Converts `x` to its inferred normalized component distribution.

        Returns
        -------
        A `x.shape[0] x n_components` tensor containing the normalized component distribution.
        """
        cell_component_unnormalized_dist = self._unnormalized_transform(x)
        return (
            cell_component_unnormalized_dist
            / cell_component_unnormalized_dist.sum(axis=1)[:, np.newaxis]
        )

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
        The negative ELBO.
        """
        return -Trace_ELBO().loss(self.model, self.guide, x, library, n_obs=n_obs)
