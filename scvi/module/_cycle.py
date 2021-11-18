# -*- coding: utf-8 -*-
"""Main module."""

import torch
import torch.nn as nn
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi import _CONSTANTS
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import Encoder, FCLayers

torch.backends.cudnn.benchmark = True


class CycleVAE(BaseModuleClass):
    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_latent_batch: int = 2,
        n_layers: int = 1,
        n_it_per_cycle_check=5,
    ):
        super().__init__()
        self.n_batch = n_batch
        self.n_latent_batch = n_latent_batch

        self.z_enc = Encoder(
            n_input,
            n_latent,
            distribution="normal",
            n_layers=n_layers,
            n_hidden=n_hidden,
        )

        self.zbatch_mean = nn.Parameter(torch.randn((n_batch, n_latent_batch)))
        self.zbatch_logvar = nn.Parameter(torch.randn((n_batch, n_latent_batch)))
        self.zbatch_activation = nn.Softplus()

        self.cycle_it = 0
        self.n_it_per_cycle_check = n_it_per_cycle_check

        self.z_full_dec = nn.Sequential(
            FCLayers(
                n_latent + n_latent_batch, n_input, n_layers=n_layers, n_hidden=n_hidden
            ),
            nn.Softmax(dim=-1),
        )

        self.px_r = nn.Parameter(torch.randn(n_input))

    def _is_cycle_check_it(self):
        return self.cycle_it % self.n_it_per_cycle_check == 0

    def _get_inference_input(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]

        input_dict = dict(x=x, batch_index=batch_index)
        return input_dict

    def _get_generative_input(self, _tensors, inference_outputs):
        z = inference_outputs["z"]
        zbatch = inference_outputs["zbatch"]
        log_library = inference_outputs["log_library"]

        input_dict = dict(
            z=z,
            zbatch=zbatch,
            log_library=log_library,
        )
        return input_dict

    @auto_move_data
    def inference(self, x, batch_index, n_samples=1, is_cycle_check=False):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = torch.log1p(x)

        log_library = torch.log(x.sum(1)).unsqueeze(1)

        qz_m, qz_v, z = self.z_enc(x_)
        batch_assignments = batch_index.long().squeeze(-1)
        qzbatch_m = self.zbatch_mean[batch_assignments]
        qzbatch_v = self.zbatch_activation(self.zbatch_logvar)[batch_assignments]
        zbatch = Normal(qzbatch_m, qzbatch_v.sqrt()).rsample()

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            z = Normal(qz_m, qz_v.sqrt()).sample((n_samples,))
            zbatch = Normal(qzbatch_m, qzbatch_v.sqrt()).sample((n_samples,))
            log_library = log_library.unsqueeze(0).expand(
                (n_samples, log_library.size(0), log_library.size(1))
            )

        return dict(
            z=z,
            qz_m=qz_m,
            qz_v=qz_v,
            zbatch=zbatch,
            qzbatch_m=qzbatch_m,
            qzbatch_v=qzbatch_v,
            log_library=log_library,
        )

    @auto_move_data
    def generative(
        self,
        z,
        zbatch,
        log_library,
    ):
        """Runs the generative model."""
        inputs = torch.cat([z, zbatch], dim=-1)
        px_scale = self.z_full_dec(inputs) * log_library.exp()

        return dict(
            px_scale=px_scale,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
    ):
        self.cycle_it += 1

        x = tensors[_CONSTANTS.X_KEY].clone()

        if self._is_cycle_check_it():
            with torch.no_grad():
                _, cycle_generative_outputs = self.forward(
                    tensors,
                    inference_kwargs=dict(is_cycle_check=True),
                    compute_loss=False,
                )
                cycle_x = NegativeBinomial(
                    cycle_generative_outputs["px_scale"], self.px_r.exp()
                ).sample()

            tensors[_CONSTANTS.X_KEY] = cycle_x
            inference_outputs, generative_outputs = self.forward(
                tensors, compute_loss=False
            )

        z = inference_outputs["z"]
        zbatch = inference_outputs["zbatch"]

        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]
        qzbatch_m = inference_outputs["qzbatch_m"]
        qzbatch_v = inference_outputs["qzbatch_v"]
        px_scale = generative_outputs["px_scale"]

        dist_px = NegativeBinomial(px_scale, self.px_r.exp())

        dist_qz = Normal(qz_m, qz_v.sqrt())
        dist_pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        dist_qzbatch = Normal(qzbatch_m, qzbatch_v.sqrt())
        dist_pzbatch = Normal(torch.zeros_like(zbatch), torch.ones_like(zbatch))

        reconstruction_loss = dist_px.log_prob(x).sum(dim=-1)
        kl_z = kl(dist_qz, dist_pz).sum(dim=-1)
        kl_zbatch = kl(dist_qzbatch, dist_pzbatch).sum(dim=-1)

        loss = (-reconstruction_loss + kl_z + kl_zbatch).mean()
        return LossRecorder(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            kl_z=kl_z,
            kl_zbatch=kl_zbatch,
        )
