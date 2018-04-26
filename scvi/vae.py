# -*- coding: utf-8 -*-
"""Main module."""
import collections

import torch
import torch.nn as nn
from torch.distributions import Normal

from scvi.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.utils import one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class VAE(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1,
                 dropout_rate=0.1, dispersion="gene", log_variational=True, kl_scale=1, reconstruction_loss="zinb",
                 batch=False, n_batch=0, using_cuda=True, n_labels=None):
        super(VAE, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_latent
        self.n_hidden = n_hidden
        self.n_input = n_input
        self.n_layers = n_layers
        self.library = 0
        self.z = 0
        self.dispersion = dispersion
        self.log_variational = log_variational
        self.kl_scale = kl_scale
        self.reconstruction_loss = reconstruction_loss
        self.n_batch = n_batch
        self.using_cuda = using_cuda and torch.cuda.is_available()
        # boolean indicating whether we want to take the batch indexes into account
        self.batch = batch

        if self.dispersion == "gene":
            self.register_buffer('px_r', torch.randn(self.n_input, ))

        self.z_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                                 dropout_rate=dropout_rate)
        self.l_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=1, n_layers=1,
                                 dropout_rate=dropout_rate)
        self.decoder = Decoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                               dropout_rate=dropout_rate, batch=batch, n_batch=n_batch)

    def sample_from_posterior_z(self, x, y=None):
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder.forward(x)
        return z

    def sample_from_posterior_l(self, x, y=None):
        # Here we compute as little as possible to have q(z|x)
        ql_m, ql_v, library = self.l_encoder.forward(x)
        return library

    def get_sample_scale(self, x, y=None, batch_index=None):
        z = self.sample_from_posterior_z(x)
        px = self.decoder.px_decoder_batch(z, batch_index)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y=None, batch_index=None):
        z = self.sample_from_posterior_z(x)
        library = self.sample_from_posterior_l(x)
        px = self.decoder.px_decoder_batch(z, batch_index)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def sample(self, z):
        return self.px_scale_decoder(z)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):  # same signature as loss
        # Parameters for z latent distribution
        if self.batch and batch_index is None:
            raise ("This VAE was trained to take batches into account:"
                   "please provide batch indexes when running the forward pass")

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


# Encoder
class Encoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1):
        super(Encoder, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_latent
        self.n_hidden = n_hidden
        self.n_input = n_input
        self.n_layers = n_layers
        # Encoding q(z/x)
        # There is always a first layer
        self.first_layer = nn.Sequential(
            nn.Dropout(p=self.dropout_rate),
            nn.Linear(n_input, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.hidden_layers = nn.Sequential(collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(
                nn.Dropout(p=self.dropout_rate),
                nn.Linear(n_hidden, n_hidden),
                nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
                nn.ReLU())) for i in range(1, n_layers)]))

        # Then, there are two different layers that compute the means and the variances of the normal distribution
        # that represents the data in the latent space
        self.encoder = nn.Sequential(self.first_layer, self.hidden_layers)
        self.mean_encoder = nn.Linear(n_hidden, n_latent)
        self.var_encoder = nn.Linear(n_hidden, n_latent)

    def reparameterize(self, mu, var):
        return Normal(mu, var.sqrt()).rsample()

    def forward(self, x):
        # Parameters for latent distribution
        q = self.encoder(x)
        q_m = self.mean_encoder(q)
        q_v = torch.exp(self.var_encoder(q))
        latent = self.reparameterize(q_m, q_v)
        return q_m, q_v, latent


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, batch=False, n_batch=0):
        super(Decoder, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_latent
        self.n_hidden = n_hidden
        self.n_input = n_input
        self.n_layers = n_layers
        self.n_batch = n_batch
        self.batch = batch

        if batch:
            self.n_hidden_real = n_hidden + n_batch
            self.n_latent_real = n_latent + n_batch
        else:
            self.n_hidden_real = n_hidden
            self.n_latent_real = n_latent

        # There is always a first layer
        self.decoder_first_layer = nn.Sequential(
            nn.Linear(self.n_latent_real, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.decoder_hidden_layers = nn.Sequential(
            collections.OrderedDict([('Layer {}'.format(i), nn.Sequential(
                nn.Dropout(p=self.dropout_rate),
                nn.Linear(self.n_hidden, n_hidden),
                nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
                nn.ReLU())) for i in range(1, n_layers)]))

        self.px_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers)

        # mean gamma
        self.px_scale_decoder = nn.Sequential(nn.Linear(self.n_hidden_real, self.n_input), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Linear(self.n_hidden_real, self.n_input)

        # dropout
        self.px_dropout_decoder = nn.Linear(self.n_hidden_real, self.n_input)

    def forward(self, dispersion, z, library, batch_index=None):
        # The decoder returns values for the parameters of the ZINB distribution
        px = self.px_decoder_batch(z, batch_index)
        px_scale = self.px_scale_decoder(px)
        px_dropout = self.px_dropout_decoder(px)
        px_rate = torch.exp(library) * px_scale
        if dispersion == "gene-cell":
            px_r = self.px_r_decoder(px)
            return px_scale, px_r, px_rate, px_dropout
        elif dispersion == "gene":
            return px_scale, px_rate, px_dropout

    def px_decoder_batch(self, z, batch_index):
        if self.batch:
            one_hot_batch = one_hot(batch_index, self.n_batch)
            z = torch.cat((z, one_hot_batch), 1)
        px = self.px_decoder(z)
        if self.batch:
            px = torch.cat((px, one_hot_batch), 1)
        return px
