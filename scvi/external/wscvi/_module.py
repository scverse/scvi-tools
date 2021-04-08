"""
This is a proof of concept to identify the simplest way to incorporate IW flavors to SCVI
TODO:
Simplify inference, generative and just work on the loss?
Simplify init as it is too complex just because of the softplus thing

So far:
in VAE
- Simplify multi sampling parts 
- Option to use Encoder with  softplus

In WVAE:
simple init
custom loss where all the log densities are computed in a separate function
custom training to be able to play with n_samples of generative/inference

Remaining interrogation:
Where/how should I code the cross evaluation part
i.e., I have x, z in some dataloader
and I want to evaluate the log ratios
"""


from typing import Iterable, Optional

import torch
from torch.distributions import Normal
import torch.nn.functional as F
import torch.nn as nn

from scvi._compat import Literal
from scvi import _CONSTANTS
from scvi.module import VAE
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import DecoderSCVI, Encoder, LinearDecoderSCVI, one_hot


def reparameterize_gaussian(mu, var):
    return Normal(mu, var.sqrt()).rsample()


class WVAE(VAE):
    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_particles: int = 1,
        loss_type: Literal["IWELBO", "ELBO"] = "IWELBO",
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: float = 0.1,
        dispersion: str = "gene",
        log_variational: bool = True,
        gene_likelihood: str = "nb",
        latent_distribution: str = "normal",
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_observed_lib_size: bool = False,
        ):
            super().__init__(
                n_input=n_input,
                n_batch=n_batch,
                n_labels=n_labels,
                n_hidden=n_hidden,
                n_latent=n_latent,
                n_layers=n_layers,
                n_continuous_cov=n_continuous_cov,
                n_cats_per_cov=n_cats_per_cov,
                dropout_rate=dropout_rate,
                dispersion=dispersion,
                log_variational=log_variational,
                gene_likelihood=gene_likelihood,
                latent_distribution=latent_distribution,
                encode_covariates=encode_covariates,
                deeply_inject_covariates=deeply_inject_covariates,
                use_batch_norm=use_batch_norm,
                use_layer_norm=use_layer_norm,
                use_observed_lib_size=use_observed_lib_size,
            )

            n_input_encoder = n_input + n_continuous_cov * encode_covariates
            cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)
            encoder_cat_list = cat_list if encode_covariates else None
            use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
            use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"

            # Only differences with VAE
            # Maybe worth to at least generalize custom encoder to lighten this part
            self.n_particles = n_particles
            self.loss_type = loss_type
            self.z_encoder = CustomEncoder(
                n_input=n_input_encoder,
                n_output=n_latent,
                n_cat_list=encoder_cat_list,
                n_layers=n_layers,
                n_hidden=n_hidden,
                dropout_rate=dropout_rate,
                distribution=latent_distribution,
                inject_covariates=deeply_inject_covariates,
                use_batch_norm=use_batch_norm_encoder,
                use_layer_norm=use_layer_norm_encoder,
            )
            # l encoder goes from n_input-dimensional data to 1-d library size
            self.l_encoder = CustomEncoder(
                n_input=n_input_encoder,
                n_output=1,
                n_layers=1,
                n_cat_list=encoder_cat_list,
                n_hidden=n_hidden,
                dropout_rate=dropout_rate,
                inject_covariates=deeply_inject_covariates,
                use_batch_norm=use_batch_norm_encoder,
                use_layer_norm=use_layer_norm_encoder,
            )

    def _get_inference_input(self, tensors):
        res =  super()._get_inference_input(tensors)
        local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
        local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
        res["local_l_mean"] = local_l_mean
        res["local_l_var"] = local_l_var
        return res

    def _get_generative_input(self, tensors, inference_outputs):
        res = super()._get_generative_input(tensors, inference_outputs)
        x = tensors[_CONSTANTS.X_KEY]
        res["x"] = x
        return res

    @auto_move_data
    def inference(
        self,
        x,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        local_l_mean=None,
        local_l_var=None,
        n_samples=None,
    ):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        if n_samples is None:
            n_samples = (self.n_particles,)
        else:
            n_samples = (n_samples,)
        x_ = x
        if self.use_observed_lib_size:
            library = torch.log(x.sum(1)).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)

        if cont_covs is not None and self.encode_covariates is True:
            encoder_input = torch.cat((x_, cont_covs), dim=-1)
        else:
            encoder_input = x_
        if cat_covs is not None and self.encode_covariates is True:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()
        # qz_m, qz_v, z = self.z_encoder(encoder_input, batch_index, *categorical_input)
        # ql_m, ql_v, library_encoded = self.l_encoder(
        #     encoder_input, batch_index, *categorical_input
        # )
        qz_m, qz_v, _ = self.z_encoder(encoder_input, batch_index, *categorical_input)
        ql_m, ql_v, _ = self.l_encoder(encoder_input, batch_index, *categorical_input)

        zdist = Normal(qz_m, qz_v.sqrt())
        untran_z = zdist.sample(n_samples)
        log_qz = zdist.log_prob(untran_z).sum(-1)
        z = self.z_encoder.z_transformation(untran_z)
        zmean = torch.zeros_like(z)
        zstd = torch.ones_like(z)
        log_pz = Normal(zmean, zstd).log_prob(z).sum(-1)

        if self.use_observed_lib_size:
            library = library.unsqueeze(0).expand(
                (n_samples, library.size(0), library.size(1))
            )
            log_ql = 0.0
            log_pl = 0.0
        else:
            ldist = Normal(ql_m, ql_v.sqrt())
            library = ldist.sample(n_samples)
            log_ql = ldist.log_prob(library).sum(-1)
            log_pl = (
                Normal(local_l_mean, torch.sqrt(local_l_var)).log_prob(library).sum(-1)
            )

        outputs = dict(
            z=z,
            qz_m=qz_m,
            qz_v=qz_v,
            ql_m=ql_m,
            ql_v=ql_v,
            library=library,
            log_ql=log_ql,
            log_qz=log_qz,
            log_pz=log_pz,
            log_pl=log_pl,
        )
        return outputs

    @auto_move_data
    def generative(
        self, z, library, batch_index, x, cont_covs=None, cat_covs=None, y=None
    ):
        """Runs the generative model."""
        decoder_input = z if cont_covs is None else torch.cat([z, cont_covs], dim=-1)
        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()
        px_scale, px_r, px_rate, px_dropout = self.decoder(
            self.dispersion, decoder_input, library, batch_index, *categorical_input, y
        )
        if self.dispersion == "gene-label":
            px_r = F.linear(
                one_hot(y, self.n_labels), self.px_r
            )  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r

        px_r = torch.exp(px_r)
        log_px_latents = -self.get_reconstruction_loss(x, px_rate, px_r, px_dropout)

        return dict(
            px_scale=px_scale,
            px_r=px_r,
            px_rate=px_rate,
            px_dropout=px_dropout,
            log_px_latents=log_px_latents,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        log_px_latents = generative_outputs["log_px_latents"]
        log_ql = inference_outputs["log_ql"]
        log_qz = inference_outputs["log_qz"]
        log_pz = inference_outputs["log_pz"]
        log_pl = inference_outputs["log_pl"]
        log_p_joint = log_px_latents + log_pz + log_pl
        log_var = log_ql + log_qz
        log_ratios = log_p_joint - log_var
        assert log_ratios.ndim == 2
        if self.loss_type == "ELBO":
            loss =  -log_ratios.mean()
        elif self.loss_type == "IWELBO":
            loss = -(torch.softmax(log_ratios, dim=0).detach() * log_ratios).sum(dim=0).mean()
        else:
            raise ValueError("Unknown loss type {}".format(self.loss_type))
        reconst_loss = log_px_latents.mean(0)
        kl_local = torch.tensor(0.)
        kl_global = torch.tensor(0.)
        return LossRecorder(loss, reconst_loss, kl_local, kl_global)

class CustomEncoder(Encoder):
    """
    Custom encoder using softplus activation for the std instead of the exponential as it showed increased stability
    """

    def __init__(self, **encoder_kwargs):
        super().__init__(**encoder_kwargs)

    def forward(self, x: torch.Tensor, *cat_list: int):
        q = self.encoder(x, *cat_list)
        q_m = self.mean_encoder(q)
        q_v = nn.Softplus()(self.var_encoder(q))
        latent = self.z_transformation(reparameterize_gaussian(q_m, q_v))
        return q_m, q_v, latent
