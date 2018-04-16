# -*- coding: utf-8 -*-
"""Main module."""
import collections

import torch
import torch.nn as nn
from torch.autograd import Variable

from scvi.log_likelihood import log_zinb_positive, log_nb_positive

torch.backends.cudnn.benchmark = True


# VAE model
class VAE(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1,
                 dropout_rate=0.1, dispersion="gene", log_variational=True, kl_scale=1, reconstruction_loss="zinb",
                 batch=False, n_batch=0):
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
        # boolean indicating whether we want to take the batch indexes into account
        self.batch = batch

        if self.dispersion == "gene":
            self.register_buffer('px_r', Variable(torch.randn(self.n_input, )))

        self.encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                               dropout_rate=dropout_rate)
        self.decoder = Decoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                               dropout_rate=dropout_rate, batch=batch, n_batch=n_batch)

    def sample_from_posterior(self, x):
        # Here we compute as little as possible to have q(z|x)
        qz = self.encoder.z_encoder(x)
        qz_m = self.encoder.z_mean_encoder(qz)
        qz_v = torch.exp(self.encoder.z_var_encoder(qz))
        return self.reparameterize(qz_m, qz_v)

    def get_sample_rate(self, x, batch_index=None):
        z = self.sample_from_posterior(x)
        px = self.decoder.px_decoder_batch(z, batch_index)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def reparameterize(self, mu, var):
        std = torch.sqrt(var)
        eps = Variable(std.data.new(std.size()).normal_())
        return eps.mul(std).add_(mu)

    def forward(self, x, batch_index=None):
        # Parameters for z latent distribution
        if self.batch and batch_index is None:
            raise ("This VAE was trained to take batches into account:"
                   "please provide batch indexes when running the forward pass")

        if self.log_variational:
            x = torch.log(1 + x)

        qz_m, qz_v, ql_m, ql_v = self.encoder.forward(x)

        # Sampling
        self.z = self.reparameterize(qz_m, qz_v)
        self.library = self.reparameterize(ql_m, ql_v)

        if self.dispersion == "gene-cell":
            px_scale, self.px_r, px_rate, px_dropout = self.decoder.forward(self.dispersion,
                                                                            self.z, self.library, batch_index)
        elif self.dispersion == "gene":
            px_scale, px_rate, px_dropout = self.decoder.forward(self.dispersion,
                                                                 self.z, self.library, batch_index)

        return px_scale, self.px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v

    def sample(self, z):
        return self.px_scale_decoder(z)

    def loss(self, sampled_batch, local_l_mean, local_l_var, kl_ponderation, batch_index=None):

        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v = self(sampled_batch, batch_index)

        # Reconstruction Loss
        if self.reconstruction_loss == 'zinb':
            reconst_loss = -log_zinb_positive(sampled_batch, px_rate, torch.exp(px_r), px_dropout)
        elif self.reconstruction_loss == 'nb':
            reconst_loss = -log_nb_positive(sampled_batch, px_rate, torch.exp(px_r))

        # KL Divergence
        kl_divergence_z = torch.sum(0.5 * (qz_m ** 2 + qz_v - torch.log(qz_v + 1e-8) - 1), dim=1)
        kl_divergence_l = torch.sum(0.5 * (
            ((ql_m - local_l_mean) ** 2) / local_l_var + ql_v / local_l_var
            + torch.log(local_l_var + 1e-8) - torch.log(ql_v + 1e-8) - 1), dim=1)

        kl_divergence = (kl_divergence_z + kl_divergence_l)

        # Train Loss # Not real total loss
        train_loss = torch.mean(reconst_loss + kl_ponderation * kl_divergence)
        return train_loss, reconst_loss, kl_divergence


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
        self.z_encoder = nn.Sequential(self.first_layer, self.hidden_layers)
        self.z_mean_encoder = nn.Linear(n_hidden, n_latent)
        self.z_var_encoder = nn.Linear(n_hidden, n_latent)

        # Encoding q(l/x)
        # The process is similar than for encoding q(z/x), except there is always only one hidden layer
        self.l_encoder = nn.Sequential(
            nn.Dropout(p=self.dropout_rate),
            nn.Linear(n_input, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        self.l_mean_encoder = nn.Linear(n_hidden, 1)
        self.l_var_encoder = nn.Linear(n_hidden, 1)

    def forward(self, x):
        # Parameters for z latent distribution
        qz = self.z_encoder(x)
        qz_m = self.z_mean_encoder(qz)
        qz_v = torch.exp(self.z_var_encoder(qz))

        # Parameters for l latent distribution
        ql = self.l_encoder(x)
        ql_m = self.l_mean_encoder(ql)
        ql_v = torch.exp(self.l_var_encoder(ql))

        return qz_m, qz_v, ql_m, ql_v


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

        def one_hot(batch_index, n_batch, dtype):
            if torch.cuda.is_available():
                batch_index = batch_index.cuda()
            onehot = batch_index.new(batch_index.size(0), n_batch).fill_(0)
            onehot.scatter_(1, batch_index, 1)
            return Variable(onehot.type(dtype))

        if self.batch:
            one_hot_batch = one_hot(batch_index, self.n_batch, z.data.type())
            z = torch.cat((z, one_hot_batch), 1)
        px = self.px_decoder(z)
        if self.batch:
            px = torch.cat((px, one_hot_batch), 1)
        return px
