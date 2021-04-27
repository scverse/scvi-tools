# -*- coding: utf-8 -*-
import numpy as np
import torch
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi import _CONSTANTS
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import Encoder, FCLayers

torch.backends.cudnn.benchmark = True


# Conditional VAE model
class VAEC(BaseModuleClass):
    """
    Conditional Variational auto-encoder model.

    This is an implementation of the CondSCVI model

    Parameters
    ----------
    n_input
        Number of input genes
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for the encoder neural network
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    """

    def __init__(
        self,
        n_input: int,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        log_variational: bool = True,
        ct_weight: np.ndarray = None,
        **module_kwargs,
    ):
        super().__init__()
        self.dispersion = "gene"
        self.n_latent = n_latent
        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.log_variational = log_variational
        self.gene_likelihood = "nb"
        self.latent_distribution = "normal"
        # Automatically deactivate if useless
        self.n_batch = 0
        self.n_labels = n_labels

        # gene dispersion
        self.px_r = torch.nn.Parameter(torch.randn(n_input))

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_cat_list=[n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=[n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
        )
        self.px_decoder = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_input), torch.nn.Softplus()
        )

        if ct_weight is not None:
            ct_weight = torch.tensor(ct_weight, dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)

    def _get_inference_input(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        y = tensors[_CONSTANTS.LABELS_KEY]

        input_dict = dict(
            x=x,
            y=y,
        )
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        y = tensors[_CONSTANTS.LABELS_KEY]

        input_dict = {
            "z": z,
            "library": library,
            "y": y,
        }
        return input_dict

    @auto_move_data
    def inference(self, x, y, n_samples=1):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        library = torch.log(x.sum(1)).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)

        qz_m, qz_v, z = self.z_encoder(x_, y)

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            # when z is normal, untran_z == z
            untran_z = Normal(qz_m, qz_v.sqrt()).sample()
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand(
                (n_samples, library.size(0), library.size(1))
            )

        outputs = dict(z=z, qz_m=qz_m, qz_v=qz_v, library=library)
        return outputs

    @auto_move_data
    def generative(self, z, library, y):
        """Runs the generative model."""
        h = self.decoder(z, y)
        px_scale = self.px_decoder(h)
        px_rate = library * px_scale

        return dict(px_scale=px_scale, px_r=self.px_r, px_rate=px_rate)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        x = tensors[_CONSTANTS.X_KEY]
        y = tensors[_CONSTANTS.LABELS_KEY]
        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]
        px_rate = generative_outputs["px_rate"]
        px_r = generative_outputs["px_r"]

        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(
            dim=1
        )

        reconst_loss = -NegativeBinomial(px_rate, logits=px_r).log_prob(x).sum(-1)
        scaling_factor = self.ct_weight[y.long()[:, 0]]
        loss = torch.mean(scaling_factor * (reconst_loss + kl_weight * kl_divergence_z))

        return LossRecorder(loss, reconst_loss, kl_divergence_z, 0.0)

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Tensors dict
        n_samples
            Number of required samples for each cell

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        inference_kwargs = dict(n_samples=n_samples)
        generative_outputs = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )[1]

        px_r = generative_outputs["px_r"]
        px_rate = generative_outputs["px_rate"]

        dist = NegativeBinomial(px_rate, logits=px_r)
        if n_samples > 1:
            exprs = dist.sample().permute(
                [1, 2, 0]
            )  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample()

        return exprs.cpu()
