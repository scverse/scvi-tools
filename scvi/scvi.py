# -*- coding: utf-8 -*-
"""Main module."""
import collections

import torch
import torch.nn as nn
from torch.autograd import Variable

from scvi.log_likelihood import log_zinb_positive

if torch.cuda.is_available():
    dtype = torch.cuda.FloatTensor
else:
    dtype = torch.FloatTensor


# VAE model
class VAE(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1,
                 dropout_rate=0.1, dispersion="gene", log_variational=True, kl_scale=1):
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
        if self.dispersion == "gene":
            self.px_r = Variable(torch.randn(self.n_input, ).type(dtype), requires_grad=False)

        self.encoder = Encoder(n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1)
        self.decoder = Decoder(n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1)

    def reparameterize(self, mu, var):
        std = torch.sqrt(var)
        eps = Variable(std.data.new(std.size()).normal_())
        return eps.mul(std).add_(mu)

    def forward(self, x):
        # Parameters for z latent distribution
        if torch.cuda.is_available():
            x = x.cuda()
        if self.log_variational:
            x = torch.log(1 + x)

        qz_m, qz_v, ql_m, ql_v = self.encoder.forward(x)

        # Sampling
        self.z = self.reparameterize(qz_m, qz_v)
        self.library = self.reparameterize(ql_m, ql_v)
        if self.dispersion == "gene-cell":
            px_scale, self.px_r, px_rate, px_dropout = self.decoder.forward(self.dispersion, self.z, self.library)
        elif self.dispersion == "gene":
            px_scale, px_rate, px_dropout = self.decoder.forward(self.dispersion, self.z, self.library)

        return px_scale, self.px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v

    def sample(self, z):
        return self.px_scale_decoder(z)

    def loss(self, sampled_batch, local_l_mean, local_l_var, kl_ponderation):
        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v = self(sampled_batch)

        # Reconstruction Loss
        reconst_loss = -log_zinb_positive(sampled_batch, px_rate, torch.exp(px_r), px_dropout)

        # KL Divergence
        kl_divergence_z = torch.sum(0.5 * (qz_m ** 2 + qz_v - torch.log(qz_v + 1e-8) - 1), dim=1)
        kl_divergence_l = torch.sum(0.5 * (
            ((ql_m - local_l_mean) ** 2) / local_l_var + ql_v / local_l_var
            + torch.log(local_l_var + 1e-8) - torch.log(ql_v + 1e-8) - 1), dim=1)

        kl_ponderation = Variable(dtype([kl_ponderation]), requires_grad=False)
        kl_divergence = (kl_divergence_z + kl_divergence_l)

        # Total Loss
        total_loss = torch.mean(reconst_loss + kl_ponderation * kl_divergence)
        return total_loss, reconst_loss, kl_divergence


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
        self.z_mean_encoder = nn.Sequential(self.first_layer, self.hidden_layers, nn.Linear(n_hidden, n_latent))
        self.z_var_encoder = nn.Sequential(self.first_layer, self.hidden_layers, nn.Linear(n_hidden, n_latent))

        # Encoding q(l/x)
        # The process is similar than for encoding q(z/x), except there is always only one hidden layer
        self.l_encoder_initial = nn.Sequential(
            nn.Dropout(p=self.dropout_rate),
            nn.Linear(n_input, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        self.l_mean_encoder = nn.Sequential(self.l_encoder_initial,
                                            nn.Linear(n_hidden, 1))
        self.l_var_encoder = nn.Sequential(self.l_encoder_initial,
                                           nn.Linear(n_hidden, 1))

        if torch.cuda.is_available():
            self.cuda()

    def forward(self, x):
        # Parameters for z latent distribution

        qz_m = self.z_mean_encoder(x)
        qz_v = torch.exp(self.z_var_encoder(x))

        # Parameters for l latent distribution
        ql_m = self.l_mean_encoder(x)
        ql_v = torch.exp(self.l_var_encoder(x))

        return qz_m, qz_v, ql_m, ql_v


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1):
        super(Decoder, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_latent
        self.n_hidden = n_hidden
        self.n_input = n_input
        self.n_layers = n_layers

        # There is always a first layer
        self.decoder_first_layer = nn.Sequential(
            nn.Linear(n_latent, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.decoder_hidden_layers = nn.Sequential(
            collections.OrderedDict([('Layer {}'.format(i), nn.Sequential(
                nn.Dropout(p=self.dropout_rate),
                nn.Linear(n_hidden, n_hidden),
                nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
                nn.ReLU())) for i in range(1, n_layers)]))

        # mean gamma
        self.px_scale_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                              nn.Linear(self.n_hidden, self.n_input), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                          nn.Linear(self.n_hidden, self.n_input))

        # dropout
        self.px_dropout_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                                nn.Linear(self.n_hidden, self.n_input))

        if torch.cuda.is_available():
            self.cuda()

    def forward(self, dispersion, z, library):
        # The decoder returns values for the parameters of the ZINB distribution
        px_scale = self.px_scale_decoder(z)
        px_dropout = self.px_dropout_decoder(z)
        px_rate = torch.exp(library) * px_scale
        if dispersion == "gene-cell":
            px_r = self.px_r_decoder(z)
            return px_scale, px_r, px_rate, px_dropout
        elif dispersion == "gene":
            return px_scale, px_rate, px_dropout
