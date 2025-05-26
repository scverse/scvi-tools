from collections.abc import Iterable

import torch
from torch import nn
from torch.distributions import Dirichlet, Normal

from scvi.nn import Decoder, FCLayers


def _identity(x):
    return x


# Encoder
class Encoder(nn.Module):
    """Encode data of ``n_input`` dimensions into a latent space of ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (data space)
    n_output
        The dimensionality of the output (latent space)
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    distribution
        Distribution of z
    var_eps
        Minimum value for the variance;
        used for numerical stability
    return_dist
        Return directly the distribution of z instead of its parameters.
    **kwargs
        Keyword args for :class:`~scvi.nn.FCLayers`
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        distribution: str = "normal",
        var_eps: float = 1e-4,
        return_dist: bool = False,
        **kwargs,
    ):
        super().__init__()

        self.distribution = distribution
        self.var_eps = var_eps
        self.encoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            **kwargs,
        )

        self.dist_encoder = nn.Linear(n_hidden, 2 * n_output)

        self.return_dist = return_dist

        if distribution == "ln":
            self.z_transformation = nn.Softmax(dim=-1)
        else:
            self.z_transformation = _identity

    def forward(self, x: torch.Tensor, *cat_list: int):
        r"""The forward computation for a single sample.

         #. Encodes the data into latent space using the encoder network
         #. Generates a mean \\( q_m \\) and variance \\( q_v \\)
         #. Samples a new value from an i.i.d. multivariate normal \\(
         # \\sim Ne(q_m, \\mathbf{I}q_v) \\)

        Parameters
        ----------
        x
            tensor with shape (n_input,)
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        3-tuple of :py:class:`torch.Tensor`
            tensors of shape ``(n_latent,)`` for mean and var, and sample

        """
        # Parameters for latent distribution
        q = self.encoder(x, *cat_list)

        q_m, q_v = self.dist_encoder(q).chunk(2, dim=-1)
        q_v = torch.nn.Softplus()(q_v) + self.var_eps

        dist = Normal(q_m, q_v.sqrt())
        latent = self.z_transformation(dist.rsample())
        if self.return_dist:
            return dist, latent
        return q_m, q_v, latent


# Decoders
class DirichletDecoder(Decoder):
    """Predict the cell type composition from the latent space

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space dimension)
    n_output
        The dimensionality of the output (number of cell type)
    n_cat_list
        A list containing the number of categories for each category of interest.
        Each category will be included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    concentration_eps
        Minimum value for the concentration; used for numerical stability
    **kwargs
        Keyword args for :class:`~scvi.nn.Decoder`
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        concentration_eps: float = 1e-6,
        **kwargs,
    ):
        super().__init__(
            n_input=n_input,
            n_output=n_output,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            **kwargs,
        )

        self.concentration_eps = concentration_eps

    def forward(self, x: torch.Tensor, *cat_list: int):
        p = self.decoder(x, *cat_list)
        p_m = self.mean_decoder(p)

        p_m = torch.nn.Softplus()(p_m) + self.concentration_eps

        dist = Dirichlet(p_m)

        return dist


class NicheDecoder(nn.Module):
    """Decodes data from latent space to niche embedding space.

    ``n_input`` dimensions to ``n_output``
    dimensions using a fully-connected neural network of ``n_hidden`` layers.
    Output is the mean and variance of a multivariate Gaussian

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space)
    n_output
        The dimensionality of the output (data space)
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    kwargs
        Keyword args for :class:`~scvi.module._base.FCLayers`
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_niche_components: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        var_eps: float = 1e-4,
        **kwargs,
    ):
        super().__init__()

        self.n_niche_components = n_niche_components
        self.n_output = n_output
        self.var_eps = var_eps

        self.decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            **kwargs,
        )

        self.dist_decoder = nn.Linear(n_hidden, 2 * n_output * n_niche_components)

    def forward(self, x: torch.Tensor, *cat_list: int):
        """The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns tensors for the mean and variance of a multivariate normal distribution

        Parameters
        ----------
        x
            tensor with shape ``(n_input,)``
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        2-tuple of :py:class:`torch.Tensor`
            Mean and variance tensors of shape ``(n_output,)``

        """
        # Parameters for latent distribution
        p = self.decoder(x, *cat_list)

        p_m, p_v = self.dist_decoder(p).chunk(2, dim=-1)
        p_v = torch.nn.Softplus()(p_v) + self.var_eps

        if p.ndim == 2:
            p_m = p_m.view(p_m.shape[0], self.n_niche_components, self.n_output)
            p_v = p_v.view(p_v.shape[0], self.n_niche_components, self.n_output)

        elif p.ndim == 3:
            p_m = p_m.view(-1, p_m.shape[1], self.n_niche_components, self.n_output)
            p_v = p_v.view(-1, p_v.shape[1], self.n_niche_components, self.n_output)

        return p_m, p_v
