from typing import Tuple

import numpy as np
import torch
from torch.distributions import NegativeBinomial, Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.compose import AbstractVAE, SCVILoss, auto_move_data


class CellAssignModule(AbstractVAE):
    """
    Model for CellAssign.

    Parameters
    ----------
    n_genes
        Number of input genes
    n_labels
        Number of input cell types
    rho
        Binary matrix of cell type markers
    **model_kwargs
        Additional kwargs
    """

    def __init__(
        self,
        n_genes: int,
        n_labels: int,
        rho: torch.Tensor,
        **model_kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_labels = n_labels

        self.register_buffer("rho", rho)

        # perform all other initialization

    def _get_inference_input(self, tensors):
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[_CONSTANTS.X_KEY]
        y = tensors[_CONSTANTS.LABELS_KEY]

        input_dict = dict(x=x, y=y)
        return input_dict

    @auto_move_data
    def inference(self):
        return {}

    @auto_move_data
    def generative(self, x, y):
        # compute phi and mean of NegBin

        # compute gamma

        # return dict(
        #     mu = mu,
        #     phi = phi,
        #     gamma = gamma,
        # )
        raise NotImplementedError

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        n_obs: int = 1.0,
    ):
        # compute Q and put it in loss
        # second term in SCVILoss is Q per cell without prior terms shape is (n_cells,)
        # third term is log prob of prior terms in Q
        # generative_outputs is a dict of the return value from `generative(...)`
        # assume that `n_obs` is the number of training data points
        # see the VAE class for how to compute NegBin log probability

        # return SCVILoss(loss, Q, prior_log_prob, 0.0)
        raise NotImplementedError

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        raise NotImplementedError("No sampling method for CellAssign")
