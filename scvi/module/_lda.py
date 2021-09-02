from typing import Dict, Iterable, Union

import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
import torch.nn.functional as F
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
        return F.softmax(self.fcalpha(h), dim=-1)


class LDAModule(PyroBaseModuleClass):
    def __init__(self, n_input: int, n_hidden: int, n_topics: int):
        super().__init__()

        self.n_input = n_input
        self.n_hidden = n_hidden
        self.n_topics = n_topics

        self.encoder = LDAEncoder(n_input, n_topics, n_hidden)

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:

        x = tensor_dict[_CONSTANTS.X_KEY]
        library = torch.sum(x, dim=1)
        return (x, library), {}

    def model(self, x: torch.Tensor, library: torch.Tensor):
        # Hyperparameters.
        with pyro.plate("topics", self.n_topics):
            alpha = pyro.sample("alpha", dist.Gamma(1.0 / self.n_topics, 1.0))
            beta = pyro.sample(
                "beta", dist.Dirichlet(torch.ones(self.n_input) / self.n_input)
            )

        # Theta generative model.
        for cell_idx in pyro.plate("cells", x.shape[0]):
            theta = pyro.sample("theta", dist.Dirichlet(alpha))

            pyro.sample(
                f"gene_counts_{cell_idx}",
                dist.Multinomial(library[cell_idx], beta @ theta),
                obs=x[cell_idx, :],
            )

    def guide(self, x: torch.Tensor, _library: torch.Tensor):
        # Hyperparameter guides.
        alpha_posterior = pyro.param(
            "alpha_posterior",
            lambda: torch.ones(self.n_topics),
            constraint=constraints.postitive,
        )
        beta_posterior = pyro.param(
            "beta_posterior",
            lambda: torch.ones(self.n_topics, self.n_input),
            constraint=constraints.greater_than(0.5),
        )
        with pyro.plate("topics", self.n_topics):
            pyro.sample("alpha", dist.Gamma(alpha_posterior, 1.0))
            pyro.sample("beta", dist.Dirichlet(beta_posterior))

        # Theta guide.
        pyro.module("encoder", self.encoder)
        with pyro.plate("cells", x.shape[0]):
            alpha = self.encoder(x)
            pyro.sample("theta", dist.Dirichlet(alpha).to_event(1))
