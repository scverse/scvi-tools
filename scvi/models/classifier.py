import torch
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from torch import nn as nn
from torch.nn import Linear

from scvi.metrics.clustering import get_latent_mean
from scvi.models.modules import FCLayers


class Classifier(nn.Module):
    def __init__(self, n_input, n_hidden=128, n_labels=10, n_layers=1, dropout_rate=0.1, use_cuda=False):
        super(Classifier, self).__init__()
        self.layers = FCLayers(n_in=n_input, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                               dropout_rate=dropout_rate)
        self.classifier = nn.Sequential(self.layers, nn.Linear(n_hidden, n_labels), nn.Softmax(dim=-1))

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.cuda()

    def forward(self, x):
        return self.classifier(x)


class LinearLogRegClassifier(Linear):
    def update_parameters(self, vae, data_loader_classification):
        data_train, _, labels_train = get_latent_mean(vae, data_loader_classification)
        log_reg = LogisticRegression()
        log_reg.fit(data_train, labels_train)

        self.weight = nn.Parameter(torch.from_numpy(log_reg.coef_).to(self.weight.device).type(self.weight.dtype))
        self.bias = nn.Parameter(torch.from_numpy(log_reg.intercept_).to(self.bias.device).type(self.bias.dtype))

    def forward(self, input):
        return F.softmax(super(LinearLogRegClassifier, self).forward(input), dim=-1)
