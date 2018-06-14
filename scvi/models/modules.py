import collections

import torch
from torch import nn as nn
from torch.distributions import Normal

from scvi.models.utils import one_hot


class FCLayers(nn.Module):
    def __init__(self, n_in, n_out, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(FCLayers, self).__init__()
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]
        self.fc_layers = FCLayers._sequential(layers_dim, dropout_rate=dropout_rate)

    def forward(self, x, *os):
        return self.fc_layers(x)

    @staticmethod
    def create(n_in, n_out, n_cat=0, n_hidden=128, n_layers=1, dropout_rate=0.1):
        if type(n_cat) is int:
            if n_cat == 0:
                return FCLayers(n_in, n_out, n_hidden=n_hidden, n_layers=n_layers, dropout_rate=dropout_rate)
            else:
                return OneHotFCLayers(n_in, n_out, n_cat=n_cat, n_hidden=n_hidden, n_layers=n_layers,
                                      dropout_rate=dropout_rate)
        elif type(n_cat) is list:
            return ManyOneHotFCLayers(n_in, n_out, n_cat_list=n_cat,
                                      n_hidden=n_hidden, n_layers=n_layers, dropout_rate=dropout_rate)

    @staticmethod
    def _sequential(layers_dim, n_cat=0, dropout_rate=0.1):
        return nn.Sequential(collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(
                nn.Dropout(p=dropout_rate),
                nn.Linear(n_in + n_cat, n_out),
                nn.BatchNorm1d(n_out, eps=1e-3, momentum=0.99),
                nn.ReLU())) for i, (n_in, n_out) in enumerate(zip(layers_dim[:-1], layers_dim[1:]))]))


class OneHotFCLayers(nn.Module):
    def __init__(self, n_in, n_out, n_cat, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(OneHotFCLayers, self).__init__()
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]
        self.n_cat = n_cat
        self.fc_layers = FCLayers._sequential(layers_dim, n_cat=n_cat, dropout_rate=dropout_rate)

    def forward(self, x, o, *os):
        if o.size(1) != self.n_cat:
            o = one_hot(o, self.n_cat)
        for layer in self.fc_layers:
            x = layer(torch.cat((x, o), 1))
        return x


class ManyOneHotFCLayers(nn.Module):
    def __init__(self, n_in, n_out, n_cat_list, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(ManyOneHotFCLayers, self).__init__()
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]
        self.n_cat_list = n_cat_list
        self.fc_layers = FCLayers._sequential(layers_dim, n_cat=sum(n_cat_list), dropout_rate=dropout_rate)

    def forward(self, x, *os):
        one_hot_os = []
        for i, o in enumerate(os):
            if o is not None and self.n_cat_list[i]:
                one_hot_o = o
                if o.size(1) != self.n_cat_list[i]:
                    one_hot_o = one_hot(o, self.n_cat_list[i])
                elif o.size(1) == 1 and self.n_cat_list[i] == 1:
                    one_hot_o = o.type(torch.float32)
                one_hot_os += [one_hot_o]
        for layer in self.fc_layers:
            x = layer(torch.cat((x,) + tuple(one_hot_os), 1))
        return x


# Encoder
class Encoder(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_latent=10, n_cat=0, n_layers=1, dropout_rate=0.1):
        super(Encoder, self).__init__()
        self.encoder = FCLayers.create(n_in=n_input, n_out=n_hidden, n_cat=n_cat, n_layers=n_layers, n_hidden=n_hidden,
                                       dropout_rate=dropout_rate)
        self.mean_encoder = nn.Linear(n_hidden, n_latent)
        self.var_encoder = nn.Linear(n_hidden, n_latent)

    def reparameterize(self, mu, var):
        return Normal(mu, var.sqrt()).rsample()

    def forward(self, x, o=None):
        # Parameters for latent distribution
        q = self.encoder(x, o)
        q_m = self.mean_encoder(q)
        q_v = torch.exp(torch.clamp(self.var_encoder(q), -5, 5))  # (computational stability safeguard)
        latent = self.reparameterize(q_m, q_v)
        return q_m, q_v, latent


# Decoder
class DecoderSCVI(nn.Module):
    def __init__(self, n_latent, n_input, n_hidden=128, n_layers=1, dropout_rate=0.1, n_batch=0, n_labels=0):
        super(DecoderSCVI, self).__init__()
        self.n_batch = n_batch
        self.px_decoder = FCLayers.create(n_in=n_latent, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                                          dropout_rate=dropout_rate, n_cat=[n_batch, n_labels])

        # mean gamma
        self.px_scale_decoder = nn.Sequential(nn.Linear(n_hidden, n_input), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Linear(n_hidden, n_input)

        # dropout
        self.px_dropout_decoder = nn.Linear(n_hidden, n_input)

    def forward(self, dispersion, z, library, batch_index=None, y=None):
        # The decoder returns values for the parameters of the ZINB distribution
        px = self.px_decoder(z, batch_index, y)
        px_scale = self.px_scale_decoder(px)
        px_dropout = self.px_dropout_decoder(px)
        # Clamp to high value: exp(12) ~ 160000 to avoid nans (computational stability)
        px_rate = torch.exp(torch.clamp(library, max=12)) * px_scale
        if dispersion == "gene-cell":
            px_r = self.px_r_decoder(px)
            return px_scale, px_r, px_rate, px_dropout
        else:  # dispersion == "gene" / "gene-batch" / "gene-label"
            return px_scale, px_rate, px_dropout


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_latent, n_output, n_cat=0, n_hidden=128, n_layers=1, dropout_rate=0.1):
        super(Decoder, self).__init__()
        self.decoder = FCLayers.create(n_in=n_latent, n_out=n_hidden, n_cat=n_cat, n_layers=n_layers,
                                       n_hidden=n_hidden, dropout_rate=dropout_rate)

        self.mean_decoder = nn.Linear(n_hidden, n_output)
        self.var_decoder = nn.Linear(n_hidden, n_output)

    def forward(self, x, o=None):
        # Parameters for latent distribution
        p = self.decoder(x, o)
        p_m = self.mean_decoder(p)
        p_v = torch.exp(self.var_decoder(p))
        return p_m, p_v
