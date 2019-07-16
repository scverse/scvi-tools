from torch import nn as nn

from scvi.models.modules import FCLayers


class Classifier(nn.Module):
    def __init__(
        self,
        n_input,
        n_hidden=128,
        n_labels=5,
        n_layers=1,
        dropout_rate=0.1,
        logits=False,
    ):
        super().__init__()
        layers = [
            FCLayers(
                n_in=n_input,
                n_out=n_hidden,
                n_layers=n_layers,
                n_hidden=n_hidden,
                dropout_rate=dropout_rate,
                use_batch_norm=True,
            ),
            nn.Linear(n_hidden, n_labels),
        ]
        if not logits:
            layers.append(nn.Softmax(dim=-1))

        self.classifier = nn.Sequential(*layers)

    def forward(self, x):
        return self.classifier(x)
