from typing import Dict, Iterable, Optional, Union

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
from pyro.nn import PyroModule
from torch.distributions import constraints

from scvi._constants import _CONSTANTS
from scvi.module.base import PyroBaseModuleClass, auto_move_data
from scvi.nn import FCLayers


class LDAPyroModel(PyroModule):
    def __init__(
        self,
        n_input: int,
        n_components: int,
        cell_component_prior: torch.Tensor,
        component_gene_prior: torch.Tensor,
    ):
        super().__init__()

        self.n_input = n_input
        self.n_components = n_components
        self.register_buffer("cell_component_prior", cell_component_prior)
        self.register_buffer("component_gene_prior", component_gene_prior)

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:

        x = tensor_dict[_CONSTANTS.X_KEY]
        library = torch.sum(x, dim=1)
        return (x, library), {}

    @auto_move_data
    def forward(self, x: torch.Tensor, library: torch.Tensor):
        # Component gene distributions.
        with pyro.plate("components", self.n_components):
            component_gene_dist = pyro.sample(
                "component_gene_dist",
                dist.Dirichlet(self.component_gene_prior),
            ).to(x.device)

        # Cell counts generation.
        max_library_size = int(torch.max(library).item())
        with pyro.plate("cells", x.shape[0]):
            cell_component_dist = pyro.sample(
                "cell_component_dist", dist.Dirichlet(self.cell_component_prior)
            )

            pyro.sample(
                "gene_counts",
                dist.Multinomial(
                    max_library_size, cell_component_dist @ component_gene_dist
                ),
                obs=x,
            )


class CellComponentDistPriorEncoder(nn.Module):
    def __init__(self, n_input: int, n_components: int, n_hidden: int):
        super().__init__()

        self.encoder = FCLayers(
            n_in=n_input,
            n_out=n_components,
            n_hidden=n_hidden,
            n_layers=1,
            inject_covariates=False,
        )

    def forward(self, x: torch.Tensor):
        return self.encoder(x).exp()


class LDAPyroGuide(PyroModule):
    def __init__(self, n_input: int, n_components: int, n_hidden: int):
        super().__init__()

        self.n_input = n_input
        self.n_components = n_components
        self.n_hidden = n_hidden

        self.encoder = CellComponentDistPriorEncoder(n_input, n_components, n_hidden)

    @auto_move_data
    def forward(self, x: torch.Tensor, _library: torch.Tensor):
        # Component gene distributions.
        component_gene_dist_posterior = pyro.param(
            "component_gene_dist_posterior",
            lambda: x.new_ones(self.n_components, self.n_input),
            constraint=constraints.greater_than(0.5),
        )
        with pyro.plate("components", self.n_components):
            pyro.sample(
                "component_gene_dist", dist.Dirichlet(component_gene_dist_posterior)
            )

        # Cell component distributions guide.
        with pyro.plate("cells", x.shape[0]):
            cell_component_prior = self.encoder(x)
            pyro.sample("cell_component_dist", dist.Dirichlet(cell_component_prior))


class LDAPyroModule(PyroBaseModuleClass):
    def __init__(
        self,
        n_input: int,
        n_components: int,
        n_hidden: int,
        cell_component_prior: Optional[np.ndarray] = None,
        component_gene_prior: Optional[np.ndarray] = None,
    ):
        super().__init__()

        self.n_input = n_input
        self.n_components = n_components
        self.n_hidden = n_hidden

        self.cell_component_prior = (
            torch.from_numpy(cell_component_prior)
            if cell_component_prior is not None
            else torch.ones(self.n_components) / self.n_components
        )
        self.component_gene_prior = (
            torch.from_numpy(component_gene_prior)
            if component_gene_prior is not None
            else torch.ones(self.n_input) / self.n_input
        )
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
