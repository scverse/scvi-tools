# -*- coding: utf-8 -*-


"""Main module."""
import collections
import torch
import torch.nn as nn
from torch.autograd import Variable


def benchmark():
    # TODO: better in terms of code logic to put in test?
    pass


# VAE model
class VAE(nn.Module):
    def __init__(self, nb_genes=100, h_dim=128, latent_dim=5, n_layers=1,
                 keep_prob=0.9):
        super(VAE, self).__init__()

        self.dropout_layer_prob = 1 - keep_prob
        self.latent_dim = latent_dim
        self.h_dim = h_dim
        self.nb_genes = nb_genes
        self.n_layers = n_layers
        self.l_m = 0
        self.z = 0
        # Encoding q(z/x)

        # There is always a first layer
        self.first_layer = nn.Sequential(
            nn.Dropout(p=self.dropout_layer_prob),
            nn.Linear(nb_genes, h_dim),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.hidden_layers = collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(nn.Dropout(p=self.dropout_layer_prob),
                                                  nn.Linear(h_dim, h_dim),
                                                  nn.ReLU()))
             for i in range(2, n_layers + 1)])

        # Then, there are two different layers that compute the means and the variances of the normal distribution
        # that represents the data in the latent space
        self.z_mean_encoder = nn.Sequential(self.first_layer, nn.Sequential(self.hidden_layers),
                                            nn.Linear(h_dim, latent_dim))
        self.z_logvar_encoder = nn.Sequential(self.first_layer, nn.Sequential(self.hidden_layers),
                                              nn.Linear(h_dim, latent_dim))

        # Encoding q(l/x)
        # The process is similar than for encoding q(z/x), except there is always only one hidden layer
        self.l_encoder_initial = nn.Sequential(
            nn.Dropout(p=self.dropout_layer_prob),
            nn.Linear(nb_genes, h_dim),
            nn.ReLU())

        self.l_mean_encoder = nn.Sequential(self.l_encoder_initial,
                                            nn.Linear(h_dim, 1))
        self.l_logvar_encoder = nn.Sequential(self.l_encoder_initial,
                                              nn.Linear(h_dim, 1))

        # Now the decoder that transforms a element of the latent spade into a potential gene-cell output

        # There is always a first layer
        self.decoder_first_layer = nn.Sequential(
            nn.Linear(latent_dim, h_dim),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.decoder_hidden_layers = nn.Sequential(
            collections.OrderedDict([('Layer {}'.format(i), nn.Sequential(nn.Dropout(p=self.dropout_layer_prob),
                                                                          nn.Linear(h_dim, h_dim),
                                                                          nn.ReLU()))
                                     for i in range(2, n_layers + 1)]))

        # mean gamma
        self.px_scale_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                              nn.Linear(self.h_dim, self.nb_genes), nn.Softmax(dim=1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                          nn.Linear(self.h_dim, self.nb_genes))

        # dropout
        self.px_dropout_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                                nn.Linear(self.h_dim, self.nb_genes))

    def reparameterize(self, mu, log_var):
        """"z = mean + eps * sigma where eps is sampled from N(0, 1)."""
        eps = Variable(torch.randn(mu.size(0), mu.size(1)))
        z = mu + eps * torch.exp(log_var / 2)  # 2 for converting variance to std
        return z

    def forward(self, x):
        # Parameters for z latent distribution
        mu_z = self.z_mean_encoder(x)
        log_var_z = self.z_logvar_encoder(x)

        # Parameters for l latent distribution
        mu_l = self.l_mean_encoder(x)
        log_var_l = self.l_logvar_encoder(x)

        # Sampling
        self.z = self.reparameterize(mu_z, log_var_z)
        self.l_m = self.reparameterize(mu_l, log_var_l)

        # The decoder returns values for the parameters of the ZINB distribution
        px_scale = self.px_scale_decoder(self.z)
        px_r = self.px_r_decoder(self.z)
        px_dropout = self.px_dropout_decoder(self.z)
        px_rate = torch.exp(self.l_m) * px_scale

        return px_scale, px_r, px_rate, px_dropout, mu_z, log_var_z, mu_l, log_var_l

    def sample(self, z):
        return self.px_scale_decoder(z)
