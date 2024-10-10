from collections.abc import Iterable, Sequence

import numpy as np
import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
import torch
import torch.nn as nn
import torch.utils.data
from torch.distributions import constraints
from torch.nn.functional import softmax, softplus

from scvi._constants import REGISTRY_KEYS
from scvi.module.base import PyroBaseModuleClass, auto_move_data

from ._components import ConditionalDenseNN


class DecipherPyroModule(PyroBaseModuleClass):
    """Decipher _decipher for single-cell data.

    Parameters
    ----------
    config : DecipherConfig or dict
        Configuration for the decipher _decipher.
    """

    def __init__(
        self,
        dim_genes: int,
        dim_v: int = 2,
        dim_z: int = 10,
        layers_v_to_z: Sequence[int] = (64,),
        layers_z_to_x: Sequence[int] = tuple(),
        beta: float = 0.1,
        prior: str = "normal",
    ):
        super().__init__()
        self.dim_v = dim_v
        self.dim_z = dim_z
        self.dim_genes = dim_genes
        self.layers_v_to_z = layers_v_to_z
        self.layers_z_to_x = layers_z_to_x
        self.beta = beta
        self.prior = prior

        self.decoder_v_to_z = ConditionalDenseNN(dim_v, layers_v_to_z, [dim_z] * 2)
        self.decoder_z_to_x = ConditionalDenseNN(dim_z, layers_z_to_x, [dim_genes])
        self.encoder_x_to_z = ConditionalDenseNN(dim_genes, [128], [dim_z] * 2)
        self.encoder_zx_to_v = ConditionalDenseNN(
            dim_genes + dim_z,
            [128],
            [dim_v, dim_v],
        )

        self.theta = None

        self._epsilon = 1e-5

        # Hack: to allow auto_move_data to infer device.
        self._dummy_param = nn.Parameter(torch.empty(0), requires_grad=False)

    @property
    def device(self):
        return self._dummy_param.device

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: dict[str, torch.Tensor]
    ) -> Iterable | dict:
        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        return (x,), {}

    @auto_move_data
    def model(self, x: torch.Tensor):
        pyro.module("decipher", self)

        self.theta = pyro.param(
            "theta",
            x.new_ones(self.dim_genes),
            constraint=constraints.positive,
        )

        with pyro.plate("batch", len(x)), poutine.scale(scale=1.0):
            with poutine.scale(scale=self.beta):
                if self.prior == "normal":
                    prior = dist.Normal(0, x.new_ones(self.dim_v)).to_event(1)
                elif self.prior == "gamma":
                    prior = dist.Gamma(0.3, x.new_ones(self.dim_v) * 0.8).to_event(1)
                else:
                    raise ValueError("Invalid prior, must be normal or gamma")
                v = pyro.sample("v", prior)

            z_loc, z_scale = self.decoder_v_to_z(v)
            z_scale = softplus(z_scale)
            z = pyro.sample("z", dist.Normal(z_loc, z_scale).to_event(1))

            mu = self.decoder_z_to_x(z)
            mu = softmax(mu, dim=-1)
            library_size = x.sum(axis=-1, keepdim=True)
            # Parametrization of Negative Binomial by the mean and inverse dispersion
            # See https://github.com/pytorch/pytorch/issues/42449
            # noinspection PyTypeChecker
            logit = torch.log(library_size * mu + self._epsilon) - torch.log(
                self.theta + self._epsilon
            )
            # noinspection PyUnresolvedReferences
            x_dist = dist.NegativeBinomial(
                total_count=self.theta + self._epsilon, logits=logit
            )
            pyro.sample("x", x_dist.to_event(1), obs=x)

    @auto_move_data
    def guide(self, x: torch.Tensor):
        pyro.module("decipher", self)
        with pyro.plate("batch", len(x)), poutine.scale(scale=1.0):
            x = torch.log1p(x)

            z_loc, z_scale = self.encoder_x_to_z(x)
            z_scale = softplus(z_scale)
            posterior_z = dist.Normal(z_loc, z_scale).to_event(1)
            z = pyro.sample("z", posterior_z)

            zx = torch.cat([z, x], dim=-1)
            v_loc, v_scale = self.encoder_zx_to_v(zx)
            v_scale = softplus(v_scale)
            with poutine.scale(scale=self.beta):
                if self.prior == "gamma":
                    posterior_v = dist.Gamma(softplus(v_loc), v_scale).to_event(1)
                elif self.prior == "normal" or self.prior == "student-normal":
                    posterior_v = dist.Normal(v_loc, v_scale).to_event(1)
                else:
                    raise ValueError("Invalid prior, must be normal or gamma")
                pyro.sample("v", posterior_v)
        return z_loc, v_loc, z_scale, v_scale

    def compute_v_z_numpy(self, x: np.array):
        """Compute decipher_v and decipher_z for a given input.

        Parameters
        ----------
        x : np.ndarray or torch.Tensor
            Input data of shape (n_cells, n_genes).

        Returns
        -------
        v : np.ndarray
            Decipher components v of shape (n_cells, dim_v).
        z : np.ndarray
            Decipher latent z of shape (n_cells, dim_z).
        """
        if type(x) == np.ndarray:
            x = torch.tensor(x, dtype=torch.float32)

        x = torch.log1p(x)
        z_loc, _ = self.encoder_x_to_z(x)
        zx = torch.cat([z_loc, x], dim=-1)
        v_loc, _ = self.encoder_zx_to_v(zx)
        return v_loc.detach().numpy(), z_loc.detach().numpy()

    def impute_gene_expression_numpy(self, x):
        if type(x) == np.ndarray:
            x = torch.tensor(x, dtype=torch.float32)
        z_loc, _, _, _ = self.guide(x)
        mu = self.decoder_z_to_x(z_loc)
        mu = softmax(mu, dim=-1)
        library_size = x.sum(axis=-1, keepdim=True)
        return (library_size * mu).detach().numpy()
