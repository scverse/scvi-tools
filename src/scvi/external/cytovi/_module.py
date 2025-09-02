from collections.abc import Callable, Iterable
from typing import Literal

import numpy as np
import torch
import torch.nn.functional as F
from torch import logsumexp, nn
from torch.distributions import Beta, Categorical, Independent, MixtureSameFamily, Normal
from torch.distributions import kl_divergence as kl

from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder, FCLayers

from ._constants import CYTOVI_REGISTRY_KEYS

torch.backends.cudnn.benchmark = True


class CytoVAE(BaseModuleClass):
    """Variational auto-encoder model for Cytometry.

    This is an implementation of the CytoVI model.

    Parameters
    ----------
    n_input
        Number of input proteins.
    n_batch
        Number of batches, if 0, no batch correction is performed. Default is 0.
    n_labels
        Number of labels. Default is 0.
    n_hidden
        Number of nodes per hidden layer. Default is 128.
    n_latent
        Dimensionality of the latent space. Default is 10.
    n_layers
        Number of hidden layers used for encoder and decoder NNs. Default is 1.
    n_continuous_cov
        Number of continuous covariates. Default is 0.
    n_cats_per_cov
        Number of categories for each extra categorical covariate. Default is None.
    dropout_rate
        Dropout rate for neural networks. Default is 0.1.
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization. Default is False.
    protein_likelihood
        One of the following protein likelihood distributions:
        * ``'normal'`` - Normal distribution
        * ``'beta'`` - Beta distribution
        Default is "normal".
    latent_distribution
        One of the following latent space distributions:
        * ``'normal'`` - Isotropic normal
        * ``'ln'`` - Logistic normal with normal params N(0, 1)
        Default is "normal".
    encode_covariates
        Whether to concatenate covariates to expression in encoder. Default is False.
    deeply_inject_covariates
        Whether to concatenate covariates into output of hidden layers in encoder/decoder.
        This option only applies when `n_layers` > 1. The covariates are concatenated to
        the input of subsequent hidden layers. Default is True.
    use_batch_norm
        Whether to use batch norm in layers. Default is "both".
    use_layer_norm
        Whether to use layer norm in layers. Default is "none".
    var_activation
        Callable used to ensure positivity of the variational distributions' variance.
        When `None`, defaults to `torch.exp`. Default is None.
    encoder_marker_mask
        List of indices to select specific markers for the encoder. Default is None.
    extra_encoder_kwargs
        Extra keyword arguments for the encoder. Default is None.
    extra_decoder_kwargs
        Extra keyword arguments for the decoder. Default is None.
    scale_activation
        Activation function for scaling factors. Default is None.
    prior_mixture
        Whether to use a mixture of gaussian prior. Default is True.
    prior_mixture_k
        Number of components in the mixture of gaussian prior. Default is 20.
    prior_label_weight
        Weight for the prior label. Default is 10.
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        dropout_rate: float = 0.1,
        log_variational: bool = False,
        protein_likelihood: Literal["normal", "beta"] = "normal",
        latent_distribution: Literal["normal", "ln"] = "normal",
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        var_activation: Callable | None = None,
        encoder_marker_mask: list | None = None,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        scale_activation: Literal["softplus", None] | None = None,
        prior_mixture: bool | None = True,
        prior_mixture_k: int = 20,
        prior_label_weight: int | None = 10,
    ):
        super().__init__()
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.protein_likelihood = protein_likelihood
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.latent_distribution = latent_distribution
        self.encode_covariates = encode_covariates
        self.encoder_marker_mask = encoder_marker_mask
        self.prior_mixture = prior_mixture
        self.prior_label_weight = prior_label_weight

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        if encoder_marker_mask is not None:
            n_input_encoder = encoder_marker_mask.sum() + n_continuous_cov * encode_covariates
        else:
            n_input_encoder = n_input + n_continuous_cov * encode_covariates

        cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)
        encoder_cat_list = cat_list if encode_covariates else None
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        self.z_encoder = Encoder(
            n_input_encoder,
            n_latent,
            n_cat_list=encoder_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=latent_distribution,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            var_activation=var_activation,
            return_dist=True,
            **_extra_encoder_kwargs,
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        n_input_decoder = n_latent + n_continuous_cov
        _extra_decoder_kwargs = extra_decoder_kwargs or {}
        self.decoder = DecoderCytoVI(
            n_input_decoder,
            n_input,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            scale_activation=scale_activation,
            protein_likelihood=protein_likelihood,
            **_extra_decoder_kwargs,
        )

        if self.prior_mixture is True:
            if self.n_labels > 1:
                prior_mixture_k = n_labels
                self.prior_means = torch.nn.Parameter(torch.zeros([prior_mixture_k, n_latent]))
            else:
                self.prior_means = torch.nn.Parameter(torch.randn([prior_mixture_k, n_latent]))
            self.prior_log_scales = torch.nn.Parameter(torch.zeros([prior_mixture_k, n_latent]))
            self.prior_logits = torch.nn.Parameter(torch.zeros([prior_mixture_k]))

    def _get_inference_input(
        self,
        tensors,
    ):
        batch_index = tensors[CYTOVI_REGISTRY_KEYS.BATCH_KEY]

        cont_key = CYTOVI_REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = CYTOVI_REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        x = tensors[CYTOVI_REGISTRY_KEYS.X_KEY]

        if self.encoder_marker_mask is not None:
            x_ = x[..., self.encoder_marker_mask]
        else:
            x_ = x
        input_dict = {
            "x": x_,
            "batch_index": batch_index,
            "cont_covs": cont_covs,
            "cat_covs": cat_covs,
        }

        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        batch_index = tensors[CYTOVI_REGISTRY_KEYS.BATCH_KEY]
        y = tensors[CYTOVI_REGISTRY_KEYS.LABELS_KEY]

        cont_key = CYTOVI_REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = CYTOVI_REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        input_dict = {
            "z": z,
            "batch_index": batch_index,
            "y": y,
            "cont_covs": cont_covs,
            "cat_covs": cat_covs,
        }
        return input_dict

    @auto_move_data
    def inference(self, x, batch_index, cont_covs=None, cat_covs=None, n_samples=1):
        """High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        if cont_covs is not None and self.encode_covariates:
            encoder_input = torch.cat((x_, cont_covs), dim=-1)
        else:
            encoder_input = x_
        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()
        qz, z = self.z_encoder(encoder_input, batch_index, *categorical_input)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
        outputs = {"z": z, "qz": qz}
        return outputs

    @auto_move_data
    def generative(
        self,
        z,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        y=None,
        transform_batch=None,
    ):
        """Runs the generative model."""
        # TODO: refactor forward function to not rely on y
        # Likelihood distribution
        if cont_covs is None:
            decoder_input = z
        elif z.dim() != cont_covs.dim():
            decoder_input = torch.cat(
                [z, cont_covs.unsqueeze(0).expand(z.size(0), -1, -1)], dim=-1
            )
        else:
            decoder_input = torch.cat([z, cont_covs], dim=-1)

        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        px_param1, px_param2 = self.decoder(
            decoder_input,
            batch_index,
            *categorical_input,
            y,
        )

        if self.protein_likelihood == "normal":
            px = Normal(loc=px_param1, scale=px_param2)
        elif self.protein_likelihood == "beta":
            px = Beta(concentration1=px_param1, concentration0=px_param2)

        if self.prior_mixture is True:
            prior_means = self.prior_means
            prior_scales = torch.exp(self.prior_log_scales) + 1e-4
            prior_logits = self.prior_logits

            if self.n_labels > 1:
                logits_labels = torch.where(
                    y < self.n_labels,
                    torch.nn.functional.one_hot(y.ravel(), num_classes=self.n_labels),
                    torch.zeros(y.shape[0], self.n_labels, device=y.device),
                )
                prior_logits = prior_logits + self.prior_label_weight * logits_labels
                prior_means = prior_means.expand(y.shape[0], -1, -1)
                prior_scales = prior_scales.expand(y.shape[0], -1, -1)

            cats = Categorical(logits=prior_logits)
            normal_dists = Independent(
                Normal(prior_means, prior_scales), reinterpreted_batch_ndims=1
            )
            pz = MixtureSameFamily(cats, normal_dists)

        else:
            pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {
            "px": px,
            "pz": pz,
        }

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Computes the loss function for the model."""
        x = tensors[CYTOVI_REGISTRY_KEYS.X_KEY]

        if CYTOVI_REGISTRY_KEYS.PROTEIN_NAN_MASK in tensors.keys():
            nan_mask = tensors[CYTOVI_REGISTRY_KEYS.PROTEIN_NAN_MASK]
        else:
            nan_mask = None

        if self.prior_mixture is True:
            z = inference_outputs["qz"].rsample(sample_shape=(10,))
            # sample x n_obs x n_latent
            kl_divergence_z = (
                inference_outputs["qz"].log_prob(z).sum(-1) - generative_outputs["pz"].log_prob(z)
            ).mean(0)

        else:
            kl_divergence_z = kl(inference_outputs["qz"], generative_outputs["pz"]).sum(dim=1)

        reconst_loss_int = -generative_outputs["px"].log_prob(x)

        # mask loss for unobserved values in batches
        if nan_mask is not None:
            reconst_loss = (reconst_loss_int * nan_mask).sum(-1)

        else:
            reconst_loss = reconst_loss_int.sum(-1)

        kl_local_for_warmup = kl_divergence_z

        weighted_kl_local = kl_weight * kl_local_for_warmup

        loss = torch.mean(reconst_loss + weighted_kl_local)

        kl_local = {
            "kl_divergence_z": kl_divergence_z,
        }
        return LossOutput(loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_local)

    @torch.inference_mode()
    def sample(
        self,
        tensors,
        n_samples=1,
    ) -> np.ndarray:
        r"""Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Tensors dict
        n_samples
            Number of required samples for each cell
        library_size
            Library size to scale samples to

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        inference_kwargs = {"n_samples": n_samples}
        (
            _,
            generative_outputs,
        ) = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )

        dist = generative_outputs["px"]

        if n_samples > 1:
            exprs = dist.sample().permute([1, 2, 0])  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample()

        return exprs.cpu()

    @torch.inference_mode()
    @auto_move_data
    def marginal_ll(self, tensors, n_mc_samples):
        """Computes the marginal log likelihood of the model."""
        sample_batch = tensors[CYTOVI_REGISTRY_KEYS.X_KEY]

        to_sum = torch.zeros(sample_batch.size()[0], n_mc_samples)

        for i in range(n_mc_samples):
            # Distribution parameters and sampled variables
            inference_outputs, _, losses = self.forward(tensors)
            qz = inference_outputs["qz"]
            z = inference_outputs["z"]

            # Reconstruction Loss
            reconst_loss = losses.dict_sum(losses.reconstruction_loss)

            # Log-probabilities
            p_z = (
                Normal(torch.zeros_like(qz.loc), torch.ones_like(qz.scale)).log_prob(z).sum(dim=-1)
            )
            p_x_zl = -reconst_loss
            q_z_x = qz.log_prob(z).sum(dim=-1)
            log_prob_sum = p_z + p_x_zl - q_z_x

            to_sum[:, i] = log_prob_sum

        batch_log_lkl = logsumexp(to_sum, dim=-1) - np.log(n_mc_samples)
        log_lkl = torch.sum(batch_log_lkl).item()
        return log_lkl


# Decoder
class DecoderCytoVI(nn.Module):
    """Decodes data from latent space of ``n_input`` dimensions into ``n_output`` dimensions.

    Uses a fully connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        Dimensionality of the input (latent space).
    n_output
        Dimensionality of the output (data space).
    n_cat_list
        List of category sizes for each categorical covariate;
        each is included via one-hot encoding.
    n_layers
        Number of fully connected hidden layers. Default: 1.
    n_hidden
        Number of nodes per hidden layer. Default: 128.
    inject_covariates
        Whether to inject covariates in each layer (``True``) or only in the
        first layer (``False``). Default: False.
    use_batch_norm
        Whether to use batch normalization in layers. Default: False.
    use_layer_norm
        Whether to use layer normalization in layers. Default: False.
    scale_activation
        Activation used for the scale/rate head (e.g., ``"softplus"`` or ``None``). Default: None.
    protein_likelihood
        Likelihood for protein expression (``"normal"`` or ``"beta"``). Default: ``"normal"``.
    decoder_param_eps
        Small epsilon added to decoder outputs for numerical stability. Default: 1e-8.

    Attributes
    ----------
    px_decoder
        Fully connected layers that decode the latent space.
    protein_likelihood
        Likelihood function used for protein expression.
    decoder_param_eps
        Small epsilon added to decoder parameters for numerical stability.
    px_param1_decoder
        Decoder module for the first output parameter.
    px_param2_decoder
        Decoder module for the second output parameter.

    Methods
    -------
    forward(z, *cat_list)
        Forward computation for a single sample.
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
        scale_activation: Literal["softplus", None] = None,
        protein_likelihood: Literal["normal", "beta"] = "normal",
        decoder_param_eps: float = 1e-6,
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
        )
        self.protein_likelihood = protein_likelihood
        self.decoder_param_eps = decoder_param_eps

        if scale_activation == "softplus":
            self.px_param1_decoder = nn.Sequential(
                nn.Linear(n_hidden, n_output),
                nn.Softplus(),
            )
        else:
            self.px_param1_decoder = nn.Linear(n_hidden, n_output)

        self.px_param2_decoder = nn.Linear(n_hidden, n_output)

    def forward(
        self,
        z: torch.Tensor,
        *cat_list: int,
    ):
        """The forward computation for a single sample.

        Decodes the data from the latent space using the decoder network.
        Returns parameters for the Normal or Beta distribution of protein expression.

        Parameters
        ----------
        z
            Tensor with shape `(n_input,)`.
        *cat_list
            List of category membership(s) for this sample.

        Returns
        -------
        Parameters for the Normal or Beta distribution of protein expression.

        """
        # The decoder returns values for the parameters of the emission distribution
        px = self.px_decoder(z, *cat_list)

        px_param1 = self.px_param1_decoder(px)
        px_param2 = self.px_param2_decoder(px)

        if self.protein_likelihood == "normal":
            px_param2 = F.softplus(px_param2) + self.decoder_param_eps

        elif self.protein_likelihood == "beta":
            px_param1 = F.softplus(px_param1) + self.decoder_param_eps
            px_param2 = F.softplus(px_param2) + self.decoder_param_eps

        return px_param1, px_param2
