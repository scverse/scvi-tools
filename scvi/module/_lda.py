from typing import Dict, Iterable, Union

import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
import torch.nn.functional as F
from pyro.nn import PyroModule
from torch.distributions import constraints

from scvi._constants import _CONSTANTS
from scvi.module.base import PyroBaseModuleClass


class LDAEncoder(nn.Module):
    def __init__(self, n_input: int, n_topics: int, n_hidden: int):
        super().__init__()

        self.fc1 = nn.Linear(n_input, n_hidden)
        self.fc2 = nn.Linear(n_hidden, n_hidden)
        self.fcalpha = nn.Linear(n_hidden, n_topics)

    def forward(self, x: torch.Tensor):
        h = F.softplus(self.fc1(x))
        h = F.softplus(self.fc2(h))
        return self.fcalpha(h).exp()


class LDAPyroModel(PyroModule):
    def __init__(self, n_input: int, n_topics: int):
        super().__init__()

        self.n_input = n_input
        self.n_topics = n_topics

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:

        x = tensor_dict[_CONSTANTS.X_KEY]
        library = torch.sum(x, dim=1)
        return (x, library), {}

    def forward(self, x: torch.Tensor, library: torch.Tensor):
        # Hyperparameters.
        with pyro.plate("topics", self.n_topics):
            alpha = pyro.sample("alpha", dist.Gamma(1.0 / self.n_topics, 1.0))
            beta = pyro.sample(
                "beta",
                dist.Dirichlet(torch.ones(self.n_input) / self.n_input),
            )

        # Full generative model.
        max_library_size = int(torch.max(library).item())
        with pyro.plate("cells", x.shape[0]):
            theta = pyro.sample("theta", dist.Dirichlet(alpha))

            pyro.sample(
                "gene_counts",
                dist.Multinomial(max_library_size, theta @ beta),
                obs=x,
            )


class LDAPyroGuide(PyroModule):
    def __init__(self, n_input: int, n_hidden: int, n_topics: int):
        super().__init__()

        self.n_input = n_input
        self.n_hidden = n_hidden
        self.n_topics = n_topics

        self.encoder = LDAEncoder(n_input, n_topics, n_hidden)

    def forward(self, x: torch.Tensor, _library: torch.Tensor):
        # Hyperparameter guides.
        alpha_posterior = pyro.param(
            "alpha_posterior",
            lambda: torch.ones(self.n_topics),
            constraint=constraints.positive,
        )
        beta_posterior = pyro.param(
            "beta_posterior",
            lambda: torch.ones(self.n_topics, self.n_input),
            constraint=constraints.greater_than(0.5),
        )
        with pyro.plate("topics", self.n_topics):
            pyro.sample("alpha", dist.Gamma(alpha_posterior, 1.0))
            pyro.sample("beta", dist.Dirichlet(beta_posterior))

        # Topic proportions guide.
        with pyro.plate("cells", x.shape[0]):
            alpha = self.encoder(x)
            pyro.sample("theta", dist.Dirichlet(alpha))


class LDAPyroModule(PyroBaseModuleClass):
    def __init__(self, n_input: int, n_hidden: int, n_topics: int):
        super().__init__()

        self.n_input = n_input
        self.n_hidden = n_hidden
        self.n_topics = n_topics

        self._model = LDAPyroModel(self.n_input, self.n_topics)
        self._guide = LDAPyroGuide(self.n_input, self.n_hidden, self.n_topics)
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide
