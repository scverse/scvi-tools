# -*- coding: utf-8 -*-


"""Main module."""
import collections

import numpy as np
import torch
import torch.nn as nn
from scipy.stats import truncnorm
from torch.autograd import Variable

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
            np.random.seed(1)
            pxr = np.random.normal(0, 1, (self.n_input,))
            self.px_r = Variable(torch.from_numpy(pxr).type(dtype))

        # Encoding q(z/x)
        # TODO: BatchNorm with params : eps=1e-3, momentum=0.99
        # (After checking that models coincide)
        # There is always a first layer
        self.first_layer = nn.Sequential(
            nn.Dropout(p=self.dropout_rate),
            nn.Linear(n_input, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.hidden_layers = collections.OrderedDict(
            [('Layer {}'.format(i), nn.Sequential(nn.Dropout(p=self.dropout_rate),
                                                  nn.Linear(n_hidden, n_hidden),
                                                  nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
                                                  nn.ReLU()))
             for i in range(2, n_layers + 1)])

        # Then, there are two different layers that compute the means and the variances of the normal distribution
        # that represents the data in the latent space
        self.z_mean_encoder = nn.Sequential(self.first_layer, nn.Sequential(self.hidden_layers),
                                            nn.Linear(n_hidden, n_latent))
        self.z_var_encoder = nn.Sequential(self.first_layer, nn.Sequential(self.hidden_layers),
                                           nn.Linear(n_hidden, n_latent))

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

        # Now the decoder that transforms a element of the latent spade into a potential gene-cell output

        # There is always a first layer
        self.decoder_first_layer = nn.Sequential(
            nn.Linear(n_latent, n_hidden),
            nn.BatchNorm1d(n_hidden, eps=1e-3, momentum=0.99),
            nn.ReLU())

        # We then add more layers if specified by the user, with a ReLU activation function
        self.decoder_hidden_layers = nn.Sequential(
            collections.OrderedDict([('Layer {}'.format(i), nn.Sequential(nn.Dropout(p=self.dropout_rate),
                                                                          nn.Linear(n_hidden, n_hidden),
                                                                          nn.BatchNorm1d(n_hidden,
                                                                                         eps=1e-3,
                                                                                         momentum=0.99),
                                                                          nn.ReLU()))
                                     for i in range(2, n_layers + 1)]))

        # TODO: remove this artificial checkpoint
        self.tmp_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers)
        # mean gamma
        self.px_scale_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                              nn.Linear(self.n_hidden, self.n_input), nn.Softmax(dim=-1))

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                          nn.Linear(self.n_hidden, self.n_input))

        # dropout
        self.px_dropout_decoder = nn.Sequential(self.decoder_first_layer, self.decoder_hidden_layers,
                                                nn.Linear(self.n_hidden, self.n_input))

        for m in self.modules():
            if isinstance(m, nn.Linear):
                # np.random.seed(1)
                # initializer=np.random.normal(0, 0.01,size=(m.in_features,m.out_features)).T
                STD = 0.01
                initializer = truncnorm.rvs(-2 * STD, 2 * STD, loc=0, scale=STD, size=(m.out_features, m.in_features))
                m.weight.data = torch.from_numpy(initializer).type(dtype)

                # m.weight.data.fill_(0.01)
        if torch.cuda.is_available():
            self.cuda()

    def reparameterize(self, mu, var):
        """"z = mean + eps * sigma where eps is sampled from N(0, 1)."""
        np.random.seed(1)
        eps = np.random.normal(0, 1, mu.shape)

        # eps = mu.data.new(mu.shape)
        # eps = Variable(eps.normal_(), requires_grad=False)
        eps = Variable(torch.from_numpy(eps).type(dtype), requires_grad=False)
        z = mu + eps * torch.sqrt(var)  # 2 for converting variance to std
        return z

    def forward(self, x):
        # Parameters for z latent distribution
        if self.log_variational:
            x = torch.log(1 + x)

        qz_m = self.z_mean_encoder(x)
        qz_v = torch.exp(self.z_var_encoder(x))

        # Parameters for l latent distribution
        ql_m = self.l_mean_encoder(x)
        ql_v = torch.exp(self.l_var_encoder(x))
        # Sampling
        self.z = self.reparameterize(qz_m, qz_v)
        self.library = self.reparameterize(ql_m, ql_v)

        # The decoder returns values for the parameters of the ZINB distribution
        px_scale = self.px_scale_decoder(self.z)
        if self.dispersion == "gene-cell":
            self.px_r = self.px_r_decoder(self.z)

        px_dropout = self.px_dropout_decoder(self.z)
        px_rate = torch.exp(self.library) * px_scale

        return px_scale, self.px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v

    def sample(self, z):
        return self.px_scale_decoder(z)
