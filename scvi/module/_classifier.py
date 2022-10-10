from torch import nn

from scvi.nn import FCLayers


class Classifier(nn.Module):
    """
    Basic fully-connected NN classifier.

    Parameters
    ----------
    n_input
        Number of input dimensions
    n_hidden
        Number of hidden nodes in hidden layer
    n_labels
        Numput of outputs dimensions
    n_layers
        Number of hidden layers
    dropout_rate
        dropout_rate for nodes
    logits
        Return logits or not
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    activation_fn
        Valid activation function from torch.nn
    """

    def __init__(
        self,
        n_input: int,
        n_hidden: int = 128,
        n_labels: int = 5,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        logits: bool = False,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        activation_fn: nn.Module = nn.ReLU,
    ):
        super().__init__()
        self.logits = logits
        layers = [
            FCLayers(
                n_in=n_input,
                n_out=n_hidden,
                n_layers=n_layers,
                n_hidden=n_hidden,
                dropout_rate=dropout_rate,
                use_batch_norm=use_batch_norm,
                use_layer_norm=use_layer_norm,
                activation_fn=activation_fn,
            ),
            nn.Linear(n_hidden, n_labels),
        ]
        if not logits:
            layers.append(nn.Softmax(dim=-1))

        self.classifier = nn.Sequential(*layers)

    def forward(self, x):
        """Forward computation."""
        return self.classifier(x)
