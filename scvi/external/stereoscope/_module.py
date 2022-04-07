from typing import Tuple

import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
import torch
from jax.nn import softplus as jaxsoftplus
from torch.distributions import NegativeBinomial
from torch.nn.functional import softplus

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data


class RNADeconv(BaseModuleClass):
    """
    Model of single-cell RNA-sequencing data for deconvolution of spatial transriptomics.

    Reimplementation of the ScModel module of Stereoscope [Andersson20]_:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_genes
        Number of input genes
    n_labels
        Number of input cell types
    **model_kwargs
        Additional kwargs
    """

    def __init__(
        self,
        n_genes: int,
        n_labels: int,
        **model_kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = n_labels

        # logit param for negative binomial
        self.px_o = torch.nn.Parameter(torch.randn(self.n_genes))
        self.W = torch.nn.Parameter(
            torch.randn(self.n_genes, self.n_labels)
        )  # n_genes, n_cell types

        if "ct_weight" in model_kwargs:
            ct_weight = torch.tensor(model_kwargs["ct_prop"], dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)

    @torch.no_grad()
    def get_params(self) -> Tuple[np.ndarray]:
        """
        Returns the parameters for feeding into the spatial data.

        Returns
        -------
        type
            list of tensor
        """
        return self.W.cpu().numpy(), self.px_o.cpu().numpy()

    def _get_inference_input(self, tensors):
        # we perform MAP here, so there is nothing to infer
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]

        input_dict = dict(x=x, y=y)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, y):
        """Simply build the negative binomial parameters for every cell in the minibatch."""
        px_scale = softplus(self.W)[:, y.long().ravel()].T  # cells per gene
        library = torch.sum(x, dim=1, keepdim=True)
        px_rate = library * px_scale
        scaling_factor = self.ct_weight[y.long().ravel()]

        return dict(
            px_scale=px_scale,
            px_o=self.px_o,
            px_rate=px_rate,
            library=library,
            scaling_factor=scaling_factor,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        x = tensors[REGISTRY_KEYS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        scaling_factor = generative_outputs["scaling_factor"]

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
        loss = torch.sum(scaling_factor * reconst_loss)

        return LossRecorder(loss, reconst_loss, torch.zeros((1,)), torch.tensor(0.0))

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for Stereoscope")


class SpatialDeconv:
    """
    Model of single-cell RNA-sequencing data for deconvolution of spatial transriptomics.

    Reimplementation of the STModel module of Stereoscope [Andersson20]_:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_spots
        Number of input spots
    sc_params
        Tuple of ndarray of shapes [(n_genes, n_labels), (n_genes)] containing the dictionnary and log dispersion parameters
    prior_weight
        Whether to sample the minibatch by the number of total observations or the monibatch size
    """

    def __init__(
        self,
        n_spots: int,
        sc_params: Tuple[np.ndarray],
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
    ):
        super().__init__()
        # unpack and copy parameters
        self.w, self.px_o = sc_params

        # setup constants
        self.n_spots = n_spots
        self.n_genes, self.n_labels = self.w.shape
        self.prior_weight = prior_weight
        self.v_init = np.random.normal(size=(self.n_labels + 1, self.n_spots))
        self.beta_init = 0.01 * np.random.normal(size=(self.n_genes, 1))

        # noise from data
        self.eta_map_init = np.random.normal(size=(self.n_genes,))

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        ind_x = tensor_dict[REGISTRY_KEYS.INDICES_KEY].astype(np.int64).squeeze()
        return (x, ind_x), {}

    def guide(self, x, ind_x):
        _, gene_plate = self.create_plates(x, ind_x)

        with gene_plate:
            numpyro.sample("eta", dist.Delta(self.eta_map_init))

    def create_plates(self, x, ind_x):
        obs_plate = numpyro.plate(
            "obs_plate", size=self.n_spots, dim=-2, subsample_size=len(ind_x)
        )
        gene_plate = numpyro.plate("gene_plate", size=self.n_genes)

        return obs_plate, gene_plate

    def model(self, x, ind_x):
        """Build the deconvolution model for every cell in the minibatch."""
        obs_plate, gene_plate = self.create_plates(x, ind_x)

        with gene_plate:
            eta = numpyro.sample("eta", dist.Normal(0, 1))

        # factor loadings
        v = numpyro.param("v", self.v_init)
        # additive gene bias
        beta = numpyro.param("beta", self.beta_init)

        beta = jaxsoftplus(beta)  # n_genes, 1
        v = jaxsoftplus(v)  # n_labels + 1, n_spots
        w = jaxsoftplus(self.w)  # n_genes, n_labels
        eps = jaxsoftplus(eta.reshape(-1, 1))  # n_genes, 1

        # account for gene specific bias and add noise
        r_hat = jnp.concatenate([beta * w, eps], axis=1)  # n_genes, n_labels + 1
        # subsample observations
        v_ind = v[:, ind_x]  # labels + 1, batch_size
        px_rate = jnp.matmul(r_hat, v_ind).transpose()  # batch_size, n_genes

        with obs_plate:
            x_dist = dist.NegativeBinomialLogits(px_rate, logits=self.px_o)
            numpyro.sample("x", x_dist.to_event(1), obs=x)
