# -*- coding: utf-8 -*-
"""Main module."""

import torch
import torch.nn as nn

from scvi.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.modules import Encoder, DecoderSCVI

torch.backends.cudnn.benchmark = True


# VAE model
class VAE(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, dispersion="gene",
                 log_variational=True, reconstruction_loss="zinb", n_batch=0, n_labels=0, using_cuda=False):
        super(VAE, self).__init__()
        self.using_cuda = using_cuda
        self.dispersion = dispersion
        self.log_variational = log_variational
        self.reconstruction_loss = reconstruction_loss
        # Automatically desactivate if useless
        self.n_batch = 0 if n_batch == 1 else n_batch

        if self.dispersion == "gene":
            self.register_buffer('px_r', torch.randn(n_input, ))

        self.z_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                                 dropout_rate=dropout_rate)
        self.l_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=1, n_layers=1,
                                 dropout_rate=dropout_rate)
        self.decoder = DecoderSCVI(n_latent, n_input, n_hidden=n_hidden, n_layers=n_layers,
                                   dropout_rate=dropout_rate, n_batch=n_batch)

    def sample_from_posterior_z(self, x, y=None):
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder.forward(x)
        return z

    def sample_from_posterior_l(self, x):
        # Here we compute as little as possible to have q(z|x)
        ql_m, ql_v, library = self.l_encoder.forward(x)
        return library

    def get_sample_scale(self, x, batch_index=None):
        z = self.sample_from_posterior_z(x)
        px = self.decoder.px_decoder_batch(z, batch_index)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y=None, batch_index=None):
        z = self.sample_from_posterior_z(x)
        library = self.sample_from_posterior_l(x)
        px = self.decoder.px_decoder(z, batch_index)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def sample(self, z):
        return self.px_scale_decoder(z)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):  # same signature as loss
        # Parameters for z latent distribution
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        # Sampling
        qz_m, qz_v, z = self.z_encoder(x_)
        ql_m, ql_v, library = self.l_encoder(x_)

        if self.dispersion == "gene-cell":
            px_scale, self.px_r, px_rate, px_dropout = self.decoder(self.dispersion, z, library, batch_index)
        elif self.dispersion == "gene":
            px_scale, px_rate, px_dropout = self.decoder(self.dispersion, z, library, batch_index)

        # Reconstruction Loss
        if self.reconstruction_loss == 'zinb':
            reconst_loss = -log_zinb_positive(x, px_rate, torch.exp(self.px_r), px_dropout)
        elif self.reconstruction_loss == 'nb':
            reconst_loss = -log_nb_positive(x, px_rate, torch.exp(self.px_r))

        # KL Divergence
        kl_divergence_z = torch.sum(0.5 * (qz_m ** 2 + qz_v - torch.log(qz_v + 1e-8) - 1), dim=1)
        kl_divergence_l = torch.sum(0.5 * (
            ((ql_m - local_l_mean) ** 2) / local_l_var + ql_v / local_l_var
            + torch.log(local_l_var + 1e-8) - torch.log(ql_v + 1e-8) - 1), dim=1)

        kl_divergence = (kl_divergence_z + kl_divergence_l)

        return reconst_loss, kl_divergence
