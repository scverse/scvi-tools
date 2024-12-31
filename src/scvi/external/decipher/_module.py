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

TEST_BETA = 1.0


class DecipherPyroModule(PyroBaseModuleClass):
    """Pyro Module for the Decipher model.

    This module implements the Decipher model for dimensionality reduction and
    interpretable representation learning in single-cell RNA sequencing data.

    Parameters
    ----------
    dim_genes
        Number of genes (features) in the dataset.
    dim_v
        Dimension of the interpretable latent space v.
    dim_z
        Dimension of the intermediate latent space z.
    layers_v_to_z
        Hidden layer sizes for the v to z decoder network.
    layers_z_to_x
        Hidden layer sizes for the z to x decoder network.
    beta
        Regularization parameter for the KL divergence.
    """

    def __init__(
        self,
        dim_genes: int,
        dim_v: int = 2,
        dim_z: int = 10,
        layers_v_to_z: Sequence[int] = (64,),
        layers_z_to_x: Sequence[int] = (),
        beta: float = 0.1,
    ):
        super().__init__()
        self.dim_v = dim_v
        self.dim_z = dim_z
        self.dim_genes = dim_genes
        self.layers_v_to_z = layers_v_to_z
        self.layers_z_to_x = layers_z_to_x
        self.beta = beta

        self.decoder_v_to_z = ConditionalDenseNN(dim_v, layers_v_to_z, [dim_z] * 2)
        self.decoder_z_to_x = ConditionalDenseNN(dim_z, layers_z_to_x, [dim_genes])
        self.encoder_x_to_z = ConditionalDenseNN(dim_genes, [128], [dim_z] * 2)
        self.encoder_zx_to_v = ConditionalDenseNN(
            dim_genes + dim_z,
            [128],
            [dim_v] * 2,
        )

        self.theta = None

        self._epsilon = 1e-5

        self.n_obs = None  # Populated by PyroTrainingPlan

        # Hack: to allow auto_move_data to infer device.
        self._dummy_param = nn.Parameter(torch.empty(0))

    @property
    def device(self):
        return self._dummy_param.device

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict: dict[str, torch.Tensor]) -> Iterable | dict:
        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        return (x,), {}

    @auto_move_data
    def model(self, x: torch.Tensor, beta: float | None = None):
        pyro.module("decipher", self)

        self.theta = pyro.param(
            "theta",
            x.new_ones(self.dim_genes),
            constraint=constraints.positive,
        )

        with (
            pyro.plate("batch", x.shape[0]),
            poutine.scale(scale=1.0),
        ):
            with poutine.scale(scale=beta or self.beta):
                prior = dist.Normal(0, x.new_ones(self.dim_v)).to_event(1)
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
            if self.device.type == "mps":
                # TODO: TORCH MPS FIX
                x_dist = dist.NegativeBinomial(
                    total_count=self.theta.contiguous() + self._epsilon, logits=logit.contiguous()
                )
            else:
                x_dist = dist.NegativeBinomial(
                    total_count=self.theta + self._epsilon, logits=logit
                )
            pyro.sample("x", x_dist.to_event(1), obs=x)

    @auto_move_data
    def guide(self, x: torch.Tensor, beta: float | None = None):
        pyro.module("decipher", self)
        with (
            pyro.plate("batch", x.shape[0]),
            poutine.scale(scale=1.0),
        ):
            x = torch.log1p(x)

            z_loc, z_scale = self.encoder_x_to_z(x)
            z_scale = softplus(z_scale)
            posterior_z = dist.Normal(z_loc, z_scale).to_event(1)
            z = pyro.sample("z", posterior_z)

            zx = torch.cat([z, x], dim=-1)
            v_loc, v_scale = self.encoder_zx_to_v(zx)
            v_scale = softplus(v_scale)
            with poutine.scale(scale=beta or self.beta):
                posterior_v = dist.Normal(v_loc, v_scale).to_event(1)
                pyro.sample("v", posterior_v)
        return z_loc, v_loc, z_scale, v_scale

    def predictive_log_likelihood(self, x: torch.Tensor, n_samples=5):
        """
        Calculate the predictive log-likelihood for a Decipher module.

        This function performs multiple runs through the dataloader to obtain
        an empirical estimate of the predictive log-likelihood. It calculates the
        log-likelihood for each run and returns the average. The beta parameter
        of the Decipher module is temporarily modified and restored even if an
        exception occurs. Used by default as an early stopping criterion.

        Parameters
        ----------
        x : torch.Tensor
            Batch of data to compute the log-likelihood for.
        n_samples : int, optional
            Number of passes through the dataloader (default is 5).

        Returns
        -------
        float
            The average estimated predictive log-likelihood across multiple runs.
        """
        log_weights = []
        for _ in range(n_samples):
            guide_trace = poutine.trace(self.guide).get_trace(x, beta=TEST_BETA)
            model_trace = poutine.trace(poutine.replay(self.model, trace=guide_trace)).get_trace(
                x, beta=TEST_BETA
            )
            log_weights.append(model_trace.log_prob_sum() - guide_trace.log_prob_sum())

        log_z = torch.logsumexp(torch.tensor(log_weights) - np.log(n_samples), 0)
        return log_z.item()
