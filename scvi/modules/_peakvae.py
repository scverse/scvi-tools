from typing import Optional, Iterable, Dict
from scvi._compat import Literal

import numpy as np
import torch
import torch.nn
from torch.distributions import Bernoulli, Normal, kl_divergence

from scvi import _CONSTANTS
from scvi.compose import (
    AbstractVAE,
    Encoder,
    FCLayers,
    SCVILoss,
    auto_move_data,
)


class Decoder(torch.nn.Module):
    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 2,
        n_hidden: int = 128,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        self.net = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            activation_fn=torch.nn.LeakyReLU,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
        )
        self.output = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_output), torch.nn.Sigmoid()
        )

    def forward(self, z: torch.Tensor, *cat_list: int):
        x = self.output(self.net(z, *cat_list))
        return x


class PEAKVAE(AbstractVAE):
    """
    Variational auto-encoder model for ATAC-seq data.

    This is an implementation of the peakVI model descibed in.

    Parameters
    ----------
    n_input_regions
        Number of input regions.
    n_batch
        Number of batches, if 0, no batch correction is performed.
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers_encoder
        Number of hidden layers used for encoder NN.
    n_layers_decoder
        Number of hidden layers used for decoder NN.
    dropout_rate
        Dropout rate for neural networks
    model_depth
        Model library size factors or not.
    region_factors
        Include region-specific factors in the model
    use_batch_norm
        One of the following

        * ``'encoder'`` - use batch normalization in the encoder only
        * ``'decoder'`` - use batch normalization in the decoder only
        * ``'none'`` - do not use batch normalization (default)
        * ``'both'`` - use batch normalization in both the encoder and decoder
    use_layer_norm
        One of the following

        * ``'encoder'`` - use layer normalization in the encoder only
        * ``'decoder'`` - use layer normalization in the decoder only
        * ``'none'`` - do not use layer normalization
        * ``'both'`` - use layer normalization in both the encoder and decoder (default)
    latent_distribution
        which latent distribution to use, options are

        * ``'normal'`` - Normal distribution (default)
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    """

    def __init__(
        self,
        n_input_regions: int,
        n_batch: int = 0,
        n_hidden: Optional[int] = None,
        n_latent: Optional[int] = None,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 2,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: float = 0.1,
        model_depth: bool = True,
        region_factors: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        latent_distribution: str = "normal",
    ):
        super().__init__()

        self.n_input_regions = n_input_regions
        self.n_hidden = int(np.sqrt(self.n_input_regions)) if n_hidden is None else n_hidden
        self.n_latent = int(np.sqrt(self.n_hidden)) if n_latent is None else n_latent
        self.n_layers_encoder = n_layers_encoder
        self.n_layers_decoder = n_layers_decoder
        self.n_cats_per_cov = n_cats_per_cov
        self.n_continuous_cov = n_continuous_cov
        self.model_depth = model_depth
        self.dropout_rate = dropout_rate
        self.latent_distribution = latent_distribution
        self.use_batch_norm_encoder = use_batch_norm in ('encoder', 'both')
        self.use_batch_norm_decoder = use_batch_norm in ('decoder', 'both')
        self.use_layer_norm_encoder = use_layer_norm in ('encoder', 'both')
        self.use_layer_norm_decoder = use_layer_norm in ('decoder', 'both')

        self.z_encoder = Encoder(
            n_input=self.n_input_regions,
            n_layers=self.n_layers_encoder,
            n_output=self.n_latent,
            n_hidden=self.n_hidden,
            dropout_rate=self.dropout_rate,
            activation_fn=torch.nn.LeakyReLU,
            distribution=self.latent_distribution, 
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
        )

        cat_list = [n_batch] + list(n_cats_per_cov) if n_cats_per_cov is not None else []
        self.z_decoder = Decoder(
            n_input=self.n_latent + self.n_continuous_cov,
            n_output=n_input_regions,
            n_hidden=self.n_hidden,
            n_cat_list=cat_list,
            n_layers=self.n_layers_decoder,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
        )

        self.d_encoder = None
        if self.model_depth:
            # Decoder class to avoid variational split
            self.d_encoder = Decoder(
                n_input=n_input_regions,
                n_output=1,
                n_hidden=self.n_hidden,
                n_layers=self.n_layers_encoder,
            )
        self.region_factors = None
        if region_factors:
            self.region_factors = torch.nn.Parameter(torch.zeros(self.n_input_regions))

    def _get_inference_input(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        input_dict = dict(
            x=x,
        )
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs, transform_batch=None):
        z = inference_outputs["z"]
        qz_m = inference_outputs["qz_m"]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]

        cont_covs = tensors.get(_CONSTANTS.CONT_COVS_KEY)

        cat_covs = tensors.get(_CONSTANTS.CAT_COVS_KEY)

        # TODO: figure out what this is
        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch
        input_dict = {
            "z": z,
            "qz_m": qz_m,
            "batch_index": batch_index,
            "cont_covs": cont_covs,
            "cat_covs": cat_covs,
        }
        return input_dict

    def get_reconstruction_loss(self, p, d, f, x):
        rl = (
            -Bernoulli(p * d * f) 
            .log_prob((x > 0).float())  # binarized data
            .sum(dim=-1)
        )
        return rl

    @auto_move_data
    def inference(
        self,
        x,
        n_samples=1,
    ) -> Dict[str, torch.Tensor]:
        """Helper function used in forward pass."""
        # Sampling
        qz_m, qz_v, z = self.z_encoder(x)
        d = self.d_encoder(x) if self.model_depth else 1

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            # when z is normal, untran_z == z
            untran_z = Normal(qz_m, qz_v.sqrt()).sample()
            z = self.z_encoder.z_transformation(untran_z)

        return dict(d=d, qz_m=qz_m, qz_v=qz_v, z=z)

    @auto_move_data
    def generative(
        self,
        z,
        qz_m,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        use_z_mean=False,
    ):
        """Runs the generative model."""

        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        latent = z if not use_z_mean else qz_m
        decoder_input = (
            latent if cont_covs is None else torch.cat([latent, cont_covs], dim=-1)
        )

        p = self.z_decoder(decoder_input, batch_index, *categorical_input)

        return dict(p=p)

    def loss(
        self, tensors, inference_outputs, generative_outputs, kl_weight: float = 1.0
    ):
        x = tensors[_CONSTANTS.X_KEY]
        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]
        d = inference_outputs["d"]
        p = generative_outputs["p"]

        kld = kl_divergence(
            Normal(qz_m, torch.sqrt(qz_v)),
            Normal(0, 1),
        ).sum(dim=1)

        f = torch.sigmoid(self.region_factors) if self.region_factors is not None else 1
        rl = self.get_reconstruction_loss(p, d, f, x)

        loss = (rl + kld * kl_weight).mean()

        return SCVILoss(loss, rl, kld, kl_global=0.0)
