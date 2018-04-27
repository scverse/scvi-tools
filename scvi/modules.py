import collections

import torch
from torch import nn as nn
from torch.distributions import Normal

from scvi.utils import one_hot


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
class DecoderSCVI(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, batch=False, n_batch=0):
        super(DecoderSCVI, self).__init__()

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


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_output=10, n_layers=1, dropout_rate=0.1):
        super(Decoder, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_output = n_output
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
        self.mean_encoder = nn.Linear(n_hidden, n_output)
        self.var_encoder = nn.Linear(n_hidden, n_output)

    def forward(self, x):
        # Parameters for latent distribution
        q = self.encoder(x)
        q_m = self.mean_encoder(q)
        q_v = torch.exp(self.var_encoder(q))
        return q_m, q_v


# Classifier
class Classifier(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_labels=10, n_layers=1, dropout_rate=0.1):
        super(Classifier, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_labels
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
        self.encoder = nn.Sequential(self.first_layer, self.hidden_layers, nn.Linear(n_hidden, n_labels),
                                     nn.Softmax(dim=-1))

    def forward(self, x):
        # Parameters for latent distribution
        return self.encoder(x)
