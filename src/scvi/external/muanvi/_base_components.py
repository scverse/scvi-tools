import torch
from torch import nn as nn
from torch.nn import functional as F

from scvi.module import Classifier
from scvi.nn import FCLayers


class Hierarchical_Classifier(nn.Module):
    """
    Hierarchical Embedding Network

    Parameters (same as Classifier )
    ----------
    n_input
        Number of input dimensions (dimensions of the latent space)
    num_classes
        number of labels in each label level in hierarchical list (ex : [2, 7])
    n_hidden
        Number of hidden nodes in one layer
    n_layers
        Number of hidden layers per NN (per independent representation)
    n_output
        Number of dimensions of each independent representation
    dropout_rate
        dropout_rate for nodes
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    concatenation
        Whether to concatenate or not the independent representations between layers
        By default no concatenation.
    """

    def __init__(
        self,
        n_input: int,
        num_classes: list,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        activation_fn: nn.Module = nn.ReLU,
        n_layers: int = 3,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
        concatenation: bool = False,
        logits: bool = True, # noqa, not used, we always return logits
    ):
        super().__init__()
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.concatenation = concatenation

        # independant representation level 1 of root level
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
            )
            for _ in range(len(num_classes))
        ]
        # neural networks to obtain independant representations of dim n_output :
        self.lvls = nn.ModuleList([nn.Sequential(layer) for layer in layers])

        self.logits = nn.ModuleList([nn.Linear(n_hidden, num_class) for num_class in num_classes])
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):
        lvl_independents = [lvl(x) for lvl in self.lvls]

        logits_level = [
            logit(lvl_independent)
            for logit, lvl_independent in zip(self.logits, lvl_independents, strict=False)
        ]
        probs_level = [self.softmax(logit_level) for logit_level in logits_level]
        return probs_level, logits_level


class HierarchicalLossNetwork(Hierarchical_Classifier):
    """
    Parameters (same as Classifier)

    ----------
    n_input
        Number of input dimensions (dimensions of the latent space)
    num_classes
        number of labels in each class in hierarchical list (ex : [2, 7])
    n_hidden
        Number of hidden nodes in one layer
    n_layers
        Number of hidden layers per NN (per independant representation)
    n_output
        Number of dimensions of each independant representation
    dropout_rate
        dropout_rate for nodes
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
        num_classes: list,
        n_hidden: int = 128,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        activation_fn: nn.Module = nn.ReLU,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        **cls_parameters,
    ):
        # initialize the Classifier
        super().__init__(
            n_input=n_input,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            activation_fn=activation_fn,
            n_layers=n_layers,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            num_classes=num_classes,
            **cls_parameters,
        )

        self.total_level = len(num_classes)

    def calculate_lloss(self, predictions, true_labels, device, weights=None):
        """
        Calculates the layer loss across all levels (multiple Cross-Entropy)

        Parameters
        ----------
        predictions
            Predictions of the model
        true labels
            Ground truth
        weights
            If we want to compute weighted cross entropy
        """
        lloss = 0
        for l in range(self.total_level):
            lloss += nn.CrossEntropyLoss()(predictions[l], true_labels[l])
        return lloss


class MultiBatchClassifier(nn.Module):
    """
    Dictionary of batch-specific fully-connected NN classifiers.

    Parameters : parameters of the batch-specific classifiers
    ----------
    n_input
        Number of input dimensions
    n_hidden
        Number of nodes in hidden layer(s). If `0`, the classifier only consists of a
        single linear layer.
    n_labels
        Numput of outputs dimensions
    n_layers
        Number of hidden layers. If `0`, the classifier only consists of a single
        linear layer.
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
    n_sites
        Number of different batch-specific classifiers to create
    **kwargs
        Keyword arguments passed into :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_sites: int,
        n_hidden: int = 128,
        n_labels: int = 5,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        activation_fn: nn.Module = nn.ReLU,
        **kwargs,
    ):
        super().__init__()
        self.n_sites = n_sites
        self.classifier_dict = nn.ModuleDict(
            {
                str(i): Classifier(
                    n_input=n_input,
                    n_hidden=n_hidden,
                    n_labels=n_labels,
                    n_layers=n_layers,
                    dropout_rate=dropout_rate,
                    use_batch_norm=use_batch_norm,
                    use_layer_norm=use_layer_norm,
                    activation_fn=activation_fn,
                    **kwargs,
                )
                for i in range(n_sites)
            }
        )

    def forward(self, x, site_index):
        """Forward computation for one mini batch of observations."""
        logits_list = []
        indices_list = []
        unique_sites = torch.unique(site_index)

        for site in unique_sites:
            site_indices = torch.nonzero(site_index == site, as_tuple=True)[0].to(x.device)
            indices_list.append(site_indices)
            x_site = x[site_indices]
            logits_list.append(self.classifier_dict[str(int(site.item()))](x_site))
        all_logits = torch.cat(logits_list, dim=0).to(x.device)
        all_indices = torch.cat(indices_list, dim=0)

        # Sort the indices to get the original minibatch order
        sorted_indices = torch.argsort(all_indices).to(x.device)
        output_logits = all_logits[sorted_indices]

        return F.softmax(output_logits, dim=-1), output_logits