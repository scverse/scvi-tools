from collections.abc import Iterable

import torch
from torch import nn
from torch.distributions import Binomial

from scvi.distributions import BetaBinomial
from scvi.external.methylvi._utils import METHYLVI_REGISTRY_KEYS
from scvi.nn import FCLayers


class BSSeqVAEMixin:
    """Shared methods for BS-seq VAE models."""

    def _compute_minibatch_reconstruction_loss(self, minibatch_size, tensors, generative_outputs):
        reconst_loss = torch.zeros(minibatch_size).to(self.device)

        for context in self.contexts:
            px_mu = generative_outputs["px_mu"][context]
            px_gamma = generative_outputs["px_gamma"][context]
            mc = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.MC_KEY}"]
            cov = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[context])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            reconst_loss += -dist.log_prob(mc).sum(dim=-1)

        return reconst_loss


class DecoderMETHYLVI(nn.Module):
    """Decodes data from latent space of ``n_input`` dimensions into ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

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
    inject_covariates
        Whether to inject covariates in each layer, or just the first (default).
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    scale_activation
        Activation layer to use for px_scale_decoder
    **kwargs
        Keyword args for :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        use_layer_norm: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.px_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        self.px_mu_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output),
            nn.Sigmoid(),
        )
        self.px_gamma_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output),
            nn.Sigmoid(),
        )

    def forward(
        self,
        dispersion: str,
        z: torch.Tensor,
        *cat_list: int,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns parameters for the beta-binomial distribution of methylation
         #. If ``dispersion != 'region-cell'`` then value for that param will be ``None``

        Parameters
        ----------
        dispersion
            One of the following

            * ``'region'`` - dispersion parameter of NB is constant per region across cells
            * ``'region-cell'`` - dispersion can differ for every region in every cell
        z :
            tensor with shape ``(n_input,)``
        library_size
            library size
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        2-tuple of :py:class:`torch.Tensor`
            parameters for the Beta distribution of mean methylation values

        """
        px = self.px_decoder(z, *cat_list)
        px_mu = self.px_mu_decoder(px)
        px_gamma = self.px_gamma_decoder(px) if dispersion == "region-cell" else None

        return px_mu, px_gamma
