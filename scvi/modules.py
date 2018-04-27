import collections

import torch
from torch import nn as nn
from torch.distributions import Normal

from scvi.utils import one_hot


class FCLayers(nn.Module):
    def __init__(self, n_in, n_out, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(FCLayers, self).__init__()

        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]
        self.fc_layers = nn.Sequential(collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(
                nn.Dropout(p=dropout_rate),
                nn.Linear(n_in, n_out),
                nn.BatchNorm1d(n_out, eps=1e-3, momentum=0.99),
                nn.ReLU())) for i, (n_in, n_out) in enumerate(zip(layers_dim[:-1], layers_dim[1:]))]))

    def forward(self, x):
        return self.fc_layers(x)


# Encoder
class Encoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1):
        super(Encoder, self).__init__()
        self.encoder = FCLayers(n_in=n_input, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                                dropout_rate=dropout_rate)
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
    def __init__(self, n_latent, n_input, n_hidden=128, n_layers=1, dropout_rate=0.1, batch=False, n_batch=0):
        super(DecoderSCVI, self).__init__()
        self.n_batch = n_batch
        self.batch = batch
        if self.batch:
            n_latent = n_latent + n_batch

        self.px_decoder = FCLayers(n_in=n_latent, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                                   dropout_rate=dropout_rate)

        # mean gamma
        self.px_scale_decoder = nn.Sequential(nn.Linear(n_hidden, n_input), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Linear(n_hidden, n_input)

        # dropout
        self.px_dropout_decoder = nn.Linear(n_hidden, n_input)

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
        return px


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_latent, n_output, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(Decoder, self).__init__()
        self.decoder = FCLayers(n_in=n_latent, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                                dropout_rate=dropout_rate)

        self.mean_decoder = nn.Linear(n_hidden, n_output)
        self.var_decoder = nn.Linear(n_hidden, n_output)

    def forward(self, x):
        # Parameters for latent distribution
        p = self.decoder(x)
        p_m = self.mean_decoder(p)
        p_v = torch.exp(self.var_decoder(p))
        return p_m, p_v


# Classifier
class Classifier(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_labels=10, n_layers=1, dropout_rate=0.1):
        super(Classifier, self).__init__()
        self.layers = FCLayers(n_in=n_input, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                               dropout_rate=dropout_rate)
        self.classifier = nn.Sequential(self.layers, nn.Linear(n_hidden, n_labels), nn.Softmax(dim=-1))

    def forward(self, x):
        # Parameters for latent distribution
        return self.classifier(x)
