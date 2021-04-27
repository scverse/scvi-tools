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
import torch.nn.functional as F
from torch.distributions import Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.module import VAE
from scvi.module.base import LossRecorder, auto_move_data
from scvi.nn import one_hot


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
        dropout_rate: float = 0.0,
        dispersion: str = "gene",
        log_variational: bool = True,
        gene_likelihood: str = "nb",
        latent_distribution: str = "normal",
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_observed_lib_size: bool = False,
        var_activation: Optional[Callable] = None,
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
            var_activation=var_activation,
        )
        self.n_particles = n_particles
        self.loss_type = loss_type

    def _get_generative_input(self, tensors, inference_outputs):
        res = super()._get_generative_input(tensors, inference_outputs)
        x = tensors[_CONSTANTS.X_KEY]
        local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
        local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
        res["local_l_mean"] = local_l_mean
        res["local_l_var"] = local_l_var
        res["x"] = x
        return res

    @auto_move_data
    def inference(
        self,
        x,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        n_samples=None,
        return_densities=True,
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
        qz_m, qz_v, _ = self.z_encoder(encoder_input, batch_index, *categorical_input)
        ql_m, ql_v, _ = self.l_encoder(encoder_input, batch_index, *categorical_input)

        zdist = Normal(qz_m, qz_v.sqrt())
        untran_z = zdist.rsample(n_samples)
        log_qz = zdist.log_prob(untran_z).sum(-1)
        z = self.z_encoder.z_transformation(untran_z)

        if self.use_observed_lib_size:
            log_ql = 0.0
            point_library = library
        else:
            ldist = Normal(ql_m, ql_v.sqrt())
            library = ldist.rsample(n_samples)
            log_ql = ldist.log_prob(library).sum(-1)
            point_library = ql_m
        log_qjoint = log_ql + log_qz
        outputs = dict(
            z=z,
            qz_m=qz_m,
            qz_v=qz_v,
            ql_m=ql_m,
            ql_v=ql_v,
            library=library,
            log_ql=log_ql,
            log_qz=log_qz,
            log_qjoint=log_qjoint,
            point_library=point_library,
        )
        return outputs

    @auto_move_data
    def generative(
        self,
        z,
        library,
        batch_index,
        x,
        local_l_mean=None,
        local_l_var=None,
        cont_covs=None,
        cat_covs=None,
        y=None,
        return_densities=True,
    ):
        """Runs the generative model."""
        assert z.ndim == 3
        decoder_input = z
        if cont_covs is not None:
            n_samples = z.shape[0]
            n_obs, n_cont = cont_covs.shape
            cont_covs_ = cont_covs.unsqueeze(0).expand(n_samples, n_obs, n_cont)
            decoder_input = torch.cat([z, cont_covs_], dim=-1)

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

        zmean = torch.zeros_like(z)
        zstd = torch.ones_like(z)
        log_pz = Normal(zmean, zstd).log_prob(z).sum(-1)
        if self.use_observed_lib_size:
            log_pl = 0.0
        else:
            log_pl = (
                Normal(local_l_mean, torch.sqrt(local_l_var)).log_prob(library).sum(-1)
            )

        log_px_latents = -self.get_reconstruction_loss(x, px_rate, px_r, px_dropout)
        log_pjoint = log_px_latents + log_pz + log_pl
        return dict(
            px_scale=px_scale,
            px_r=px_r,
            px_rate=px_rate,
            px_dropout=px_dropout,
            log_pl=log_pl,
            log_px_latents=log_px_latents,
            log_pz=log_pz,
            log_pjoint=log_pjoint,
        )

    def generative_evaluate(self, tensors, inference_outputs):
        gen_ins = self._get_generative_input(
            tensors=tensors, inference_outputs=inference_outputs
        )
        return self.generative(**gen_ins)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        log_ratios = generative_outputs["log_pjoint"] - inference_outputs["log_qjoint"]
        assert log_ratios.ndim == 2
        if self.loss_type == "ELBO":
            loss = -log_ratios.mean()
        elif self.loss_type == "IWELBO":
            loss = (
                -(torch.softmax(log_ratios, dim=0).detach() * log_ratios)
                .sum(dim=0)
                .mean()
            )
        else:
            raise ValueError("Unknown loss type {}".format(self.loss_type))
        reconst_loss = -generative_outputs["log_px_latents"].mean(0)
        kl_local = torch.tensor(0.0)
        kl_global = torch.tensor(0.0)
        return LossRecorder(loss, reconst_loss, kl_local, kl_global)
