import collections

import torch
from torch import nn as nn
from torch.distributions import Normal

from scvi.models.utils import one_hot


class FCLayers(nn.Module):
    def __init__(self, n_in, n_out, n_cat_list=[], n_layers=1, n_hidden=128, dropout_rate=0.1):
        super(FCLayers, self).__init__()
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]
        self.n_cat_list = [n_cat if n_cat > 1 else 0 for n_cat in n_cat_list]  # n_cat = 1 will be ignored
        self.fc_layers = nn.Sequential(collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(
                nn.Linear(n_in + sum(self.n_cat_list), n_out),
                nn.BatchNorm1d(n_out, eps=1e-3, momentum=0.99),
                nn.ReLU(),
                nn.Dropout(p=dropout_rate))) for i, (n_in, n_out) in enumerate(zip(layers_dim[:-1], layers_dim[1:]))]))

    def forward(self, x, *cat_list):
        one_hot_cat_list = []  # for generality in this list many indices useless.
        assert len(self.n_cat_list) <= len(cat_list), "nb. categorical args provided doesn't match init. params."
        for n_cat, cat in zip(self.n_cat_list, cat_list):
            assert not (n_cat and cat is None), "cat not provided while n_cat != 0 in init. params."
            if n_cat > 1:  # n_cat = 1 will be ignored - no additional information
                if cat.size(1) != n_cat:
                    one_hot_cat = one_hot(cat, n_cat)
                else:
                    one_hot_cat = cat  # cat has already been one_hot encoded
                one_hot_cat_list += [one_hot_cat]
        for layers in self.fc_layers:
            for layer in layers:
                if isinstance(layer, nn.Linear):
                    x = torch.cat((x, *one_hot_cat_list), 1)
                x = layer(x)
        return x


# Encoder
class Encoder(nn.Module):
    def __init__(self, n_input, n_output, n_cat_list=[], n_layers=1, n_hidden=128, dropout_rate=0.1):
        super(Encoder, self).__init__()
        self.encoder = FCLayers(n_in=n_input, n_out=n_hidden, n_cat_list=n_cat_list, n_layers=n_layers,
                                n_hidden=n_hidden, dropout_rate=dropout_rate)
        self.mean_encoder = nn.Linear(n_hidden, n_output)
        self.var_encoder = nn.Linear(n_hidden, n_output)

    def reparameterize(self, mu, var):
        return Normal(mu, var.sqrt()).rsample()

    def forward(self, x, *cat_list):
        # Parameters for latent distribution
        q = self.encoder(x, *cat_list)
        q_m = self.mean_encoder(q)
        q_v = torch.exp(torch.clamp(self.var_encoder(q), -5, 5))  # (computational stability safeguard)
        latent = self.reparameterize(q_m, q_v)
        return q_m, q_v, latent


# Decoder
class DecoderSCVI(nn.Module):
    def __init__(self, n_input, n_output, n_cat_list=[], n_layers=1, n_hidden=128, dropout_rate=0.1):
        super(DecoderSCVI, self).__init__()
        self.px_decoder = FCLayers(n_in=n_input, n_out=n_hidden, n_cat_list=n_cat_list, n_layers=n_layers,
                                   n_hidden=n_hidden, dropout_rate=dropout_rate)

        # mean gamma
        self.px_scale_decoder = nn.Sequential(nn.Linear(n_hidden, n_output), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Linear(n_hidden, n_output)

        # dropout
        self.px_dropout_decoder = nn.Linear(n_hidden, n_output)

    def forward(self, dispersion, z, library, *cat_list):
        # The decoder returns values for the parameters of the ZINB distribution
        px = self.px_decoder(z, *cat_list)
        px_scale = self.px_scale_decoder(px)
        px_dropout = self.px_dropout_decoder(px)
        # Clamp to high value: exp(12) ~ 160000 to avoid nans (computational stability)
        px_rate = torch.exp(torch.clamp(library, max=12)) * px_scale
        px_r = self.px_r_decoder(px) if dispersion == "gene-cell" else None
        return px_scale, px_r, px_rate, px_dropout


# Decoder
class Decoder(nn.Module):
    def __init__(self, n_input, n_output, n_cat_list=[], n_layers=1, n_hidden=128, dropout_rate=0.1):
        super(Decoder, self).__init__()
        self.decoder = FCLayers(n_in=n_input, n_out=n_hidden, n_cat_list=n_cat_list, n_layers=n_layers,
                                n_hidden=n_hidden, dropout_rate=dropout_rate)

        self.mean_decoder = nn.Linear(n_hidden, n_output)
        self.var_decoder = nn.Linear(n_hidden, n_output)

    def forward(self, x, *cat_list):
        # Parameters for latent distribution
        p = self.decoder(x, *cat_list)
        p_m = self.mean_decoder(p)
        p_v = torch.exp(self.var_decoder(p))
        return p_m, p_v
