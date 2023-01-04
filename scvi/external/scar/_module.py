from typing import Literal

import torch
from torch import nn
from torch.distributions import Binomial, Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial, Poisson, ZeroInflatedNegativeBinomial
from scvi.module._vae import VAE
from scvi.module.base import LossOutput, auto_move_data
from scvi.nn import FCLayers

torch.backends.cudnn.benchmark = True


class tanh(nn.Module):
    """Hyperbolic tangent activation function."""

    def __init__(self):
        super().__init__()

    def forward(self, input_x):
        """Forward pass."""
        var_tanh = torch.tanh(input_x)
        output = (1 + var_tanh) / 2
        return output


class hnormalization(nn.Module):
    """Hyperbolic normalization."""

    def __init__(self):
        super().__init__()

    def forward(self, input_x):
        """Forward pass."""
        return input_x / (input_x.sum(dim=-1, keepdim=True) + 1e-5)


class softplus(nn.Module):
    """Softplus activation function."""

    def __init__(self, sparsity=0.9):
        super().__init__()
        self.sparsity = sparsity

    def forward(self, input_x):
        """Forward pass."""
        return self._softplus(input_x)

    def _softplus(self, input_x):
        """Customized softplus activation, output range: [0, inf)"""
        var_sp = nn.functional.softplus(input_x)
        threshold = nn.functional.softplus(
            torch.tensor(-(1 - self.sparsity) * 10.0, device=input_x.device)
        )
        var_sp = var_sp - threshold
        zero = torch.zeros_like(threshold)
        var_out = torch.where(var_sp <= zero, zero, var_sp)
        return var_out


class DecoderSCAR(nn.Module):
    """
    Decodes data from latent space of ``n_input`` dimensions into ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space)
    n_output
        The dimensionality of the output (data space)
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    scale_activation
        Activation layer to use for px_scale_decoder
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_layers: int = 2,
        n_hidden: int = 150,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        scale_activation: Literal["softmax", "softplus", "softplus_sp"] = "softplus_sp",
        sparsity: float = 0.9,
    ):
        super().__init__()
        self.px_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
        )

        # mean gamma
        if scale_activation == "softmax":
            px_scale_activation = nn.Softmax(dim=-1)
        elif scale_activation == "softplus":
            px_scale_activation = nn.Softplus()
        elif scale_activation == "softplus_sp":
            px_scale_activation = softplus(sparsity)
        self.px_scale_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output),
            px_scale_activation,
            hnormalization(),
        )

        # noise ratio
        self.px_noise_decoder = nn.Sequential(
            nn.Linear(n_hidden, 1),
            tanh(),
        )

        # dropout
        self.px_dropout_decoder = nn.Linear(n_hidden, n_output)

    def forward(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
    ):
        """
        The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns parameters for the ZINB distribution of expression

        Parameters
        ----------
        z :
            tensor with shape ``(n_input,)``
        library_size
            library size

        Returns
        -------
        4-tuple of :py:class:`torch.Tensor`
            parameters for the ZINB distribution of native expression and noise ratio

        """
        # The decoder returns values for the parameters of the ZINB distribution
        px = self.px_decoder(z)
        px_scale = self.px_scale_decoder(px)
        px_dropout = self.px_dropout_decoder(px)

        # noise ratio
        px_noise_ratio = self.px_noise_decoder(px)

        # Clamp to high value: exp(12) ~ 160000 to avoid nans (computational stability)
        px_rate = torch.exp(library) * px_scale  # torch.clamp( , max=12)

        return px_scale, px_noise_ratio, px_rate, px_dropout


class SCAR_VAE(VAE):
    """
    Slightly modified version of scVI's VAE model to enable ambient RNA removal in scRNA-seq data.

    Parameters
    ----------
    ambient_profile
        The probability of occurrence of each ambient transcript.
    n_input
        Number of input genes
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for neural networks
    sparsity
        The sparsity of expected native signals. It varies between datasets,
        e.g. if one prefilters genes -- use only highly variable genes --
        the sparsity should be low; on the other hand, it should be set high
        in the case of unflitered genes.
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood
        One of
        * ``'b'`` - Binomial distribution
        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of
        * ``'normal'`` - Isotropic normal
        * ``'ln'`` - Logistic normal with normal params N(0, 1)
    use_layer_norm
        Whether to use layer norm in layers
    use_size_factor_key
        Use size_factor AnnDataField defined by the user as scaling factor in mean of conditional distribution.
        Takes priority over `use_observed_lib_size`.
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution
    library_log_means
        1 x n_batch array of means of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    library_log_vars
        1 x n_batch array of variances of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    var_activation
        Callable used to ensure positivity of the variational distributions' variance.
        When `None`, defaults to `torch.exp`.
    """

    def __init__(
        self,
        ambient_profile: torch.tensor,
        n_input: int,
        n_hidden: int = 150,
        n_latent: int = 15,
        n_layers: int = 2,
        dropout_rate: float = 0.0,
        scale_activation: Literal["softmax", "softplus", "softplus_sp"] = "softplus_sp",
        sparsity: float = 0.9,
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb", "b", "poisson"] = "b",
        latent_distribution: str = "normal",
        use_observed_lib_size: bool = True,
        **vae_kwargs,
    ):
        super().__init__(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            use_observed_lib_size=use_observed_lib_size,
            **vae_kwargs,
        )
        self.sparsity = sparsity
        self.ambient_profile = ambient_profile

        # decoder goes from n_latent-dimensional space to n_input-d data
        self.decoder = DecoderSCAR(
            n_latent,
            n_input,
            n_layers=n_layers,
            n_hidden=n_hidden,
            scale_activation=scale_activation,
            sparsity=self.sparsity,
        )

    @auto_move_data
    def generative(
        self,
        z,
        library,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        size_factor=None,
        y=None,
        transform_batch=None,
    ):
        """Runs the generative model."""
        # Likelihood distribution
        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        if not self.use_size_factor_key:
            size_factor = library

        px_scale, px_noise_ratio, px_rate, px_dropout = self.decoder(
            z,
            size_factor,
        )
        px_r = self.px_r

        # remove estimated noise
        px_scale = px_scale * (1 - px_noise_ratio)
        px_rate = px_rate * (1 - px_noise_ratio)
        pamb_scale = self.ambient_profile.to(px_scale.device) * px_noise_ratio
        pamb_rate = pamb_scale * torch.exp(size_factor)
        px_r = torch.exp(px_r)

        if self.gene_likelihood == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px_rate,
                theta=px_r,
                zi_logits=px_dropout,
                scale=px_scale,
            )
        elif self.gene_likelihood == "nb":
            px = NegativeBinomial(mu=px_rate, theta=px_r, scale=px_scale)
        elif self.gene_likelihood == "b":
            px = Binomial(total_count=torch.exp(size_factor).int(), probs=px_scale)
        elif self.gene_likelihood == "poisson":
            px = Poisson(rate=px_rate, scale=px_scale)

        # Priors
        if self.use_observed_lib_size:
            pl = None
        else:
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)
            pl = Normal(local_library_log_means, local_library_log_vars.sqrt())
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))
        return dict(
            px=px,
            pl=pl,
            pz=pz,
            pamb_scale=pamb_scale,
            pamb_rate=pamb_rate,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Compute the loss function for the model."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        kl_divergence_z = kl(inference_outputs["qz"], generative_outputs["pz"]).sum(
            dim=1
        )
        if not self.use_observed_lib_size:
            kl_divergence_l = kl(
                inference_outputs["ql"],
                generative_outputs["pl"],
            ).sum(dim=1)
        else:
            kl_divergence_l = 0.0

        # need to add the ambient rate and scale to the distribution for the loss
        px = generative_outputs["px"]
        if self.gene_likelihood == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px.mu + generative_outputs["pamb_rate"],
                theta=px.theta,
                zi_logits=px.zi_logits,
                scale=px.scale + generative_outputs["pamb_scale"],
            )
        elif self.gene_likelihood == "nb":
            px = NegativeBinomial(
                mu=px.mu + generative_outputs["pamb_rate"],
                theta=px.theta,
                scale=px.scale + generative_outputs["pamb_scale"],
            )
        elif self.gene_likelihood == "b":
            px = Binomial(
                total_count=px.total_count,
                probs=px.probs + generative_outputs["pamb_scale"],
            )
        elif self.gene_likelihood == "poisson":
            px = Poisson(
                rate=px.rate + generative_outputs["pamb_rate"],
                scale=px.scale + generative_outputs["pamb_scale"],
            )
        reconst_loss = -px.log_prob(x).sum(-1)

        kl_local_for_warmup = kl_divergence_z
        kl_local_no_warmup = kl_divergence_l

        weighted_kl_local = kl_weight * kl_local_for_warmup + kl_local_no_warmup

        loss = torch.mean(reconst_loss + weighted_kl_local)

        kl_local = dict(
            kl_divergence_l=kl_divergence_l, kl_divergence_z=kl_divergence_z
        )
        return LossOutput(
            loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_local
        )
