import collections

import torch
import torch.nn as nn
from torch.autograd import Variable

from scvi.utils import enumerate_discrete, one_hot
from scvi.vae import VAE


# VAE model - for classification: VAEC
class VAEC(VAE):
    def __init__(self, n_input, n_labels=7, y_prior=None, **kargs):
        super(VAEC, self).__init__(n_input=n_input + n_labels, **kargs)
        self.n_labels = n_labels
        self.n_input = n_input

        self.y_prior = y_prior if y_prior is not None else (1 / self.n_labels) * torch.ones(self.n_labels)
        self.classifier = Classifier(self.n_input, self.n_hidden, self.n_labels, self.n_layers, self.dropout_rate,
                                     self.using_cuda)

    def classify(self, x):
        return self.classifier(x)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        # Prepare for sampling
        xs, ys = (x, y)

        # Enumerate choices of label
        if not is_labelled:
            ys = enumerate_discrete(xs, self.n_labels)
            xs = xs.repeat(self.n_labels, 1)
            if batch_index is not None:
                batch_index = batch_index.repeat(self.n_labels, 1)
            local_l_var = local_l_var.repeat(self.n_labels, 1)
            local_l_mean = local_l_mean.repeat(self.n_labels, 1)
        else:
            ys = one_hot(ys, self.n_labels, xs.data.type())

        reconst_loss, kl_divergence = super(VAEC, self).forward(torch.cat((xs, ys), 1), local_l_mean, local_l_var,
                                                                batch_index=batch_index)

        if is_labelled:
            return reconst_loss, kl_divergence

        reconst_loss = reconst_loss.view(self.n_labels, -1)
        kl_divergence = kl_divergence.view(self.n_labels, -1)
        proba = self.classify(x)
        reconst_loss = (reconst_loss.t() * proba).sum(dim=1)
        kl_divergence = (kl_divergence.t() * proba).sum(dim=1)
        y_prior = Variable(self.y_prior.type(proba.data.type()))
        kl_divergence += torch.sum(torch.mul(proba, torch.log(y_prior) - torch.log(proba + 1e-8)), dim=-1)

        return reconst_loss, kl_divergence


# Classifier
class Classifier(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_labels=10, n_layers=1, dropout_rate=0.1, using_cuda=True):
        super(Classifier, self).__init__()

        self.dropout_rate = dropout_rate
        self.n_latent = n_labels
        self.n_hidden = n_hidden
        self.n_input = n_input
        self.n_layers = n_layers
        self.using_cuda = using_cuda and torch.cuda.is_available()
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
