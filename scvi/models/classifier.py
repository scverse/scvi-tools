from torch import nn as nn

from scvi.models.modules import FCLayers
from scvi.models.utils import Module


class Classifier(Module):
    def __init__(self, n_input, n_hidden=128, n_labels=10, n_layers=1, dropout_rate=0.1):
        super(Classifier, self).__init__()
        self.layers = FCLayers(n_in=n_input, n_out=n_hidden, n_layers=n_layers, n_hidden=n_hidden,
                               dropout_rate=dropout_rate)
        self.classifier = nn.Sequential(self.layers, nn.Linear(n_hidden, n_labels), nn.Softmax(dim=-1))

    def forward(self, x):
        return self.classifier(x)
