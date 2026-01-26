from collections.abc import Iterable

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal

from scvi.nn import FCLayers
from scvi.utils import dependencies

EPS = 1e-8


class DecoderSinglePathway(nn.Module):
    """Single-pathway decoder for unimodal likelihoods.

    Outputs a single set of parameters (scale, r, rate, dropout).
    Use for: nb, zinb, normal, lognormal, log1pnormal, ziln, gamma, zig.

    Parameters
    ----------
    n_output
        Number of output features.
    n_batches
        Number of batches.
    normalize
        If True, apply softmax normalization to output scale (for count data).
        If False, output raw scale values (for continuous data).
        Default is True for backward compatibility.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
        normalize: bool = True,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches
        self.normalize = normalize

        self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output))
        self.bias = nn.Parameter(torch.zeros(n_batches, n_output))
        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output))

        self.px_dropout_param = nn.Parameter(torch.randn(n_output) * 0.01)

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass through the decoder.

        Parameters
        ----------
        u
            Cell latent representation, shape (n_cells, n_latent).
        l
            Library size (log-scale), shape (n_cells, 1).
        batch_index
            Batch indices for cells, shape (n_cells,) or (n_cells, 1).
        v
            Feature latent representation, shape (n_features, n_latent).

        Returns
        -------
        tuple of 4 :class:`torch.Tensor`
            - px_scale: Normalized scale parameters, shape (n_cells, n_features).
            - px_r: Overdispersion parameter, shape (n_cells, n_features).
            - px_rate: Rate parameters (mean), shape (n_cells, n_features).
            - px_dropout: Dropout/zero-inflation probabilities, shape (n_features,).
        """
        # Check batch index dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias.shape[0] - 1}]"
            )

        # Get batch-specific parameters
        scale = F.softplus(self.scale_lin[batch_index])
        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        # Compute decoder output parameters
        raw_px_scale = scale * (u @ v.T) + bias

        if self.normalize:
            # Count data: softmax produces proportions, library scales to counts
            px_scale = torch.softmax(raw_px_scale, dim=-1)
            px_rate = torch.exp(l) * px_scale
        else:
            # Continuous data: no normalization, scale is positive, rate is unbounded
            px_scale = F.softplus(raw_px_scale) + EPS
            px_rate = raw_px_scale

        px_dropout = F.softplus(self.px_dropout_param)
        px_r = log_theta

        return px_scale, px_r, px_rate, px_dropout


class DecoderDualPathway(nn.Module):
    """Dual-pathway decoder for mixture likelihoods.

    Outputs two sets of parameters for foreground/background mixture models.
    Use for: nbmixture.

    Parameters
    ----------
    n_output
        Number of output features.
    n_batches
        Number of batches.
    normalize
        If True, apply softmax normalization to output scales (for count data).
        If False, output raw scale values (for continuous data).
        Default is True for backward compatibility.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
        normalize: bool = True,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches
        self.normalize = normalize

        self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output))

        self.bias1 = nn.Parameter(torch.zeros(n_batches, n_output))
        self.bias2 = nn.Parameter(torch.zeros(n_batches, n_output))

        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output))

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor,
    ) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor, tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
        """Forward pass through the decoder.

        Parameters
        ----------
        u
            Cell latent representation, shape (n_cells, n_latent).
        l
            Library size (log-scale), shape (n_cells, 1).
        batch_index
            Batch indices for cells, shape (n_cells,) or (n_cells, 1).
        v
            Feature latent representation, shape (n_features, n_latent).

        Returns
        -------
        tuple of 4 elements
            - scales: Tuple of (foreground, background) scale parameters.
              Each has shape (n_cells, n_features).
            - px_r: Overdispersion parameter, shape (n_cells, n_features).
            - rates: Tuple of (foreground, background) rate parameters (mean).
              Each has shape (n_cells, n_features).
            - mixture_logits: Logits for mixture component probabilities, shape (n_cells, n_features).
        """
        # Check batch index dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias1.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias1.shape[0] - 1}]"
            )

        # Get batch-specific parameters
        scale = F.softplus(self.scale_lin[batch_index])
        bias1 = self.bias1[batch_index]
        bias2 = self.bias2[batch_index]
        log_theta = self.log_theta[batch_index]

        # Compute decoder output parameters for both pathways
        raw_px_scale_1 = scale * (u @ v.T) + bias1
        raw_px_scale_2 = scale * (u @ v.T) + bias2

        if self.normalize:
            # Count data: softmax produces proportions, library scales to counts
            px_scale_1 = torch.softmax(raw_px_scale_1, dim=-1)
            px_scale_2 = torch.softmax(raw_px_scale_2, dim=-1)
            px_rate_1 = torch.exp(l) * px_scale_1
            px_rate_2 = torch.exp(l) * px_scale_2
        else:
            # Continuous data: no normalization, scales are positive, rates are unbounded
            px_scale_1 = F.softplus(raw_px_scale_1) + EPS
            px_scale_2 = F.softplus(raw_px_scale_2) + EPS
            px_rate_1 = raw_px_scale_1
            px_rate_2 = raw_px_scale_2

        mixture_logits = raw_px_scale_1 - raw_px_scale_2
        px_r = log_theta

        return (px_scale_1, px_scale_2), px_r, (px_rate_1, px_rate_2), mixture_logits


class DecoderProtein(nn.Module):
    """Protein decoder with optional background/foreground separation.

    Decodes latent representations to protein expression parameters using
    neural networks for flexible foreground scaling and mixture modeling.

    Parameters
    ----------
    n_input
        Dimensionality of the input (cell latent space).
    n_output_proteins
        Number of protein features.
    n_batches
        Number of batches.
    n_cat_list
        List of categorical covariate dimensions.
    dropout_rate
        Dropout rate for hidden layers.
    use_batch_norm
        Whether to use batch normalization.
    use_layer_norm
        Whether to use layer normalization.
    n_hidden
        Number of hidden units in FC layers.
    n_layers
        Number of hidden layers.
    common_scale
        If True, use shared scale parameters for background/foreground.
        If False, use separate parameters for each.
    """
    
    def __init__(
        self,
        n_input: int,
        n_output_proteins: int,
        n_batches: int,
        n_cat_list: Iterable[int] = None,
        dropout_rate: float = 0,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        n_hidden: int = 256,
        n_layers: int = 1,
        common_scale: bool = True,
    ):
        super().__init__()
        self.n_output_proteins = n_output_proteins
        self.n_batches = n_batches
        self.common_scale = common_scale

        if common_scale:
            self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output_proteins))
            self.bias = nn.Parameter(torch.zeros(n_batches, n_output_proteins))

        else:
            self.scale_lin_back = nn.Parameter(torch.zeros(n_batches, n_output_proteins))
            self.bias_back = nn.Parameter(torch.zeros(n_batches, n_output_proteins))

            self.scale_lin_fore = nn.Parameter(torch.zeros(n_batches, n_output_proteins))
            self.bias_fore = nn.Parameter(torch.zeros(n_batches, n_output_proteins))

        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output_proteins))

        linear_args = {
            "n_layers": 1,
            "use_activation": False,
            "use_batch_norm": False,
            "use_layer_norm": False,
            "dropout_rate": 0,
        }

        # Foreground scale network
        self.py_fore_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
        )
        # Predict increment factor making fg > bg
        self.py_fore_scale_decoder = FCLayers(
            n_in=n_hidden + n_input,
            n_out=n_output_proteins,
            n_cat_list=n_cat_list,
            n_layers=1,
            use_activation=True,
            use_batch_norm=False,
            use_layer_norm=False,
            dropout_rate=0,
            activation_fn=nn.ReLU,
        )

        # Mixing probability network
        self.sigmoid_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
        )
        # Predict mixing probability - (1 - protein_mixing) = fg probability
        self.py_background_decoder = FCLayers(
            n_in=n_hidden + n_input,
            n_out=n_output_proteins,
            n_cat_list=n_cat_list,
            **linear_args,
        )

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor
    ) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor, tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
        """Forward pass through the protein decoder.

        Parameters
        ----------
        u
            Cell latent representation, shape (n_cells, n_latent).
        l
            Library size (log-scale), shape (n_cells, 1).
        batch_index
            Batch indices for cells, shape (n_cells,) or (n_cells, 1).
        v
            Feature latent representation, shape (n_proteins, n_latent).

        Returns
        -------
        tuple of 4 elements
            - scales: Tuple of (background, foreground) scale parameters.
              Each has shape (n_cells, n_proteins).
            - px_r: Overdispersion parameter, shape (n_cells, n_proteins).
            - rates: Tuple of (background, foreground) rate parameters.
              Each has shape (n_cells, n_proteins).
            - mixing: Logits for background/foreground mixture weights, shape (n_cells, n_proteins).
        """
        # Check batch index dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        py_ = {}
        py_["r"] = self.log_theta[batch_index]  # px_r

        if self.common_scale:
            scale = F.softplus(self.scale_lin[batch_index])
            bias = self.bias[batch_index]

            # Parametrize the background mean using the feature/cell matrix product
            raw_px_scale = scale * (u @ v.T) + bias
            py_["scale_back"] = torch.softmax(raw_px_scale, dim=-1)
            # for fg different act function (positive + 1)
            py_["rate_back"] = torch.exp(l) * py_["scale_back"]  # calculate mean

            # Learn foreground scaling factor with a NN
            py_fore = self.py_fore_decoder(u, batch_index)
            py_fore_cat_z = torch.cat([py_fore, u], dim=-1)
            py_["scale_fore"] = self.py_fore_scale_decoder(py_fore_cat_z, batch_index) + 1 + 1e-8
            py_["rate_fore"] = py_["rate_back"] * py_["scale_fore"]

        else:
            scale_back = F.softplus(self.scale_lin_back[batch_index])
            bias_back = self.bias_back[batch_index]

            scale_fore = F.softplus(self.scale_lin_fore[batch_index])
            bias_fore = self.bias_fore[batch_index]

            # Parametrize the background mean using the feature/cell matrix product
            raw_px_scale = scale_back * (u @ v.T) + bias_back
            py_["scale_back"] = torch.softmax(raw_px_scale, dim=-1)
            # for fg different act function (positive + 1)
            py_["rate_back"] = torch.exp(l) * py_["scale_back"]  # calculate mean

            # Parametrize the background mean using the feature/cell matrix product
            raw_px_scale = scale_fore * (u @ v.T) + bias_fore
            activation_func = nn.ReLU()
            py_["scale_fore"] = activation_func(raw_px_scale) + 1 + 1e-8
            # For fg different act function (positive + 1)
            py_["rate_fore"] = torch.exp(l) * py_["scale_fore"]  # calculate mean

        # Learn the mixing logits with a NN
        p_mixing = self.sigmoid_decoder(u, batch_index)
        p_mixing_cat_z = torch.cat([p_mixing, u], dim=-1)

        py_["mixing"] = self.py_background_decoder(p_mixing_cat_z, batch_index)

        return (
            (py_["scale_back"], py_["scale_fore"]),
            py_["r"],
            (py_["rate_back"], py_["rate_fore"]),
            py_["mixing"],
        )


class GraphEncoder_glue(nn.Module):
    """Graph convolutional encoder for graph-structured data.

    Uses a Graph Convolutional Network (GCN) to encode node features and
    outputs a latent distribution.

    Parameters
    ----------
    vnum
        Number of nodes (vertices) in the graph.
    out_features
        Number of output features / latent dimension.
    """

    @dependencies("torch_geometric")
    def __init__(self, vnum: int, out_features: int):
        import torch_geometric

        super().__init__()
        self.vrepr = nn.Parameter(torch.zeros(vnum, out_features))
        self.conv = torch_geometric.nn.GCNConv(out_features, out_features)
        self.loc = nn.Linear(out_features, out_features)
        self.std_lin = nn.Linear(out_features, out_features)

    def forward(self, edge_index: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass through the graph encoder.

        Parameters
        ----------
        edge_index
            Edge indices for the graph, shape (2, n_edges).

        Returns
        -------
        tuple of 3 :class:`torch.Tensor`
            - z: Sampled latent representation from the normal distribution, shape (n_nodes, out_features).
            - mu: Mean of the latent distribution, shape (n_nodes, out_features).
            - logvar: Log variance of the latent distribution, shape (n_nodes, out_features).
        """
        h = self.conv(self.vrepr, edge_index)
        loc = self.loc(h)
        std = F.softplus(self.std_lin(h)) + EPS

        dist = Normal(loc, std)
        mu = dist.loc
        std = dist.scale
        logvar = torch.log(std**2)
        z = dist.rsample()

        return z, mu, logvar
