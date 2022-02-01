from typing import Dict, Iterable, Optional

import numpy as np
import torch
from torch.distributions import Normal, Poisson
from torch.distributions import kl_divergence as kld

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.module._peakvae import Decoder as DecoderPeakVI
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import DecoderSCVI, Encoder, FCLayers


class LibrarySizeEncoder(torch.nn.Module):
    def __init__(
        self,
        n_input: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 2,
        n_hidden: int = 128,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
        deep_inject_covariates: bool = False,
    ):
        super().__init__()
        self.px_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            activation_fn=torch.nn.LeakyReLU,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            inject_covariates=deep_inject_covariates,
        )
        self.output = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, 1), torch.nn.LeakyReLU()
        )

    def forward(self, x: torch.Tensor, *cat_list: int):
        return self.output(self.px_decoder(x, *cat_list))


class MULTIVAE(BaseModuleClass):
    """
    Variational auto-encoder model for joint paired + unpaired RNA-seq and ATAC-seq data.

    Parameters
    ----------
    n_input_regions
        Number of input regions.
    n_input_genes
        Number of input genes.
    n_batch
        Number of batches, if 0, no batch correction is performed.
    n_labels
        Number of labels, if 0, all cells are assumed to have the same label
    gene_likelihood
        The distribution to use for gene expression data. One of the following
        * ``'zinb'`` - Zero-Inflated Negative Binomial
        * ``'nb'`` - Negative Binomial
        * ``'poisson'`` - Poisson
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers_encoder
        Number of hidden layers used for encoder NN.
    n_layers_decoder
        Number of hidden layers used for decoder NN.
    dropout_rate
        Dropout rate for neural networks
    region_factors
        Include region-specific factors in the model
    use_batch_norm
        One of the following
        * ``'encoder'`` - use batch normalization in the encoder only
        * ``'decoder'`` - use batch normalization in the decoder only
        * ``'none'`` - do not use batch normalization
        * ``'both'`` - use batch normalization in both the encoder and decoder
    use_layer_norm
        One of the following
        * ``'encoder'`` - use layer normalization in the encoder only
        * ``'decoder'`` - use layer normalization in the decoder only
        * ``'none'`` - do not use layer normalization
        * ``'both'`` - use layer normalization in both the encoder and decoder
    latent_distribution
        which latent distribution to use, options are
        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    deeply_inject_covariates
        Whether to deeply inject covariates into all layers of the decoder. If False,
        covariates will only be included in the input layer.
    encode_covariates
        If True, include covariates in the input to the encoder.
    """

    ## TODO: replace n_input_regions and n_input_genes with a gene/region mask (we don't dictate which comes forst or that they're even contiguous)
    def __init__(
        self,
        n_input_regions: int = 0,
        n_input_genes: int = 0,
        n_batch: int = 0,
        n_labels: int = 0,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        n_hidden: Optional[int] = None,
        n_latent: Optional[int] = None,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 2,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: float = 0.1,
        region_factors: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        latent_distribution: str = "normal",
        deeply_inject_covariates: bool = False,
        encode_covariates: bool = False,
    ):
        super().__init__()

        # INIT PARAMS
        self.n_input_regions = n_input_regions
        self.n_input_genes = n_input_genes
        self.n_hidden = (
            int(np.sqrt(self.n_input_regions + self.n_input_genes))
            if n_hidden is None
            else n_hidden
        )
        self.n_batch = n_batch
        self.n_labels = n_labels

        self.gene_likelihood = gene_likelihood
        self.latent_distribution = latent_distribution

        self.n_latent = int(np.sqrt(self.n_hidden)) if n_latent is None else n_latent
        self.n_layers_encoder = n_layers_encoder
        self.n_layers_decoder = n_layers_decoder
        self.n_cats_per_cov = n_cats_per_cov
        self.n_continuous_cov = n_continuous_cov
        self.dropout_rate = dropout_rate

        self.use_batch_norm_encoder = use_batch_norm in ("encoder", "both")
        self.use_batch_norm_decoder = use_batch_norm in ("decoder", "both")
        self.use_layer_norm_encoder = use_layer_norm in ("encoder", "both")
        self.use_layer_norm_decoder = use_layer_norm in ("decoder", "both")
        self.encode_covariates = encode_covariates
        self.deeply_inject_covariates = deeply_inject_covariates

        cat_list = (
            [n_batch] + list(n_cats_per_cov) if n_cats_per_cov is not None else []
        )

        n_input_encoder_acc = (
            self.n_input_regions + n_continuous_cov * encode_covariates
        )
        n_input_encoder_exp = self.n_input_genes + n_continuous_cov * encode_covariates
        encoder_cat_list = cat_list if encode_covariates else None

        ## accessibility encoder
        self.z_encoder_accessibility = Encoder(
            n_input=n_input_encoder_acc,
            n_layers=self.n_layers_encoder,
            n_output=self.n_latent,
            n_hidden=self.n_hidden,
            n_cat_list=encoder_cat_list,
            dropout_rate=self.dropout_rate,
            activation_fn=torch.nn.LeakyReLU,
            distribution=self.latent_distribution,
            var_eps=0,
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
        )

        ## expression encoder
        self.z_encoder_expression = Encoder(
            n_input=n_input_encoder_exp,
            n_layers=self.n_layers_encoder,
            n_output=self.n_latent,
            n_hidden=self.n_hidden,
            n_cat_list=encoder_cat_list,
            dropout_rate=self.dropout_rate,
            activation_fn=torch.nn.LeakyReLU,
            distribution=self.latent_distribution,
            var_eps=0,
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
        )

        # expression decoder
        self.z_decoder_expression = DecoderSCVI(
            self.n_latent + self.n_continuous_cov,
            n_input_genes,
            n_cat_list=cat_list,
            n_layers=n_layers_decoder,
            n_hidden=self.n_hidden,
            inject_covariates=self.deeply_inject_covariates,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
        )

        # accessibility decoder
        self.z_decoder_accessibility = DecoderPeakVI(
            n_input=self.n_latent + self.n_continuous_cov,
            n_output=n_input_regions,
            n_hidden=self.n_hidden,
            n_cat_list=cat_list,
            n_layers=self.n_layers_decoder,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )

        ## accessibility region-specific factors
        self.region_factors = None
        if region_factors:
            self.region_factors = torch.nn.Parameter(torch.zeros(self.n_input_regions))

        ## expression dispersion parameters
        self.px_r = torch.nn.Parameter(torch.randn(n_input_genes))

        ## expression library size encoder
        self.l_encoder_expression = LibrarySizeEncoder(
            n_input_encoder_exp,
            n_cat_list=encoder_cat_list,
            n_layers=self.n_layers_encoder,
            n_hidden=self.n_hidden,
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )

        ## accessibility library size encoder
        self.l_encoder_accessibility = DecoderPeakVI(
            n_input=n_input_encoder_acc,
            n_output=1,
            n_hidden=self.n_hidden,
            n_cat_list=encoder_cat_list,
            n_layers=self.n_layers_encoder,
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )

    def _get_inference_input(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        cont_covs = tensors.get(_CONSTANTS.CONT_COVS_KEY)
        cat_covs = tensors.get(_CONSTANTS.CAT_COVS_KEY)
        input_dict = dict(
            x=x,
            batch_index=batch_index,
            cont_covs=cont_covs,
            cat_covs=cat_covs,
        )
        return input_dict

    @auto_move_data
    def inference(
        self,
        x,
        batch_index,
        cont_covs,
        cat_covs,
        n_samples=1,
    ) -> Dict[str, torch.Tensor]:

        # Get Data and Additional Covs
        x_rna = x[:, : self.n_input_genes]
        x_chr = x[:, self.n_input_genes :]

        mask_expr = x_rna.sum(dim=1) > 0
        mask_acc = x_chr.sum(dim=1) > 0

        if cont_covs is not None and self.encode_covariates:
            encoder_input_expression = torch.cat((x_rna, cont_covs), dim=-1)
            encoder_input_accessibility = torch.cat((x_chr, cont_covs), dim=-1)
        else:
            encoder_input_expression = x_rna
            encoder_input_accessibility = x_chr

        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        # Z Encoders
        qzm_acc, qzv_acc, z_acc = self.z_encoder_accessibility(
            encoder_input_accessibility, batch_index, *categorical_input
        )
        qzm_expr, qzv_expr, z_expr = self.z_encoder_expression(
            encoder_input_expression, batch_index, *categorical_input
        )

        # L encoders
        libsize_expr = self.l_encoder_expression(
            encoder_input_expression, batch_index, *categorical_input
        )
        libsize_acc = self.l_encoder_accessibility(
            encoder_input_accessibility, batch_index, *categorical_input
        )

        # ReFormat Outputs
        if n_samples > 1:
            qzm_acc = qzm_acc.unsqueeze(0).expand(
                (n_samples, qzm_acc.size(0), qzm_acc.size(1))
            )
            qzv_acc = qzv_acc.unsqueeze(0).expand(
                (n_samples, qzv_acc.size(0), qzv_acc.size(1))
            )
            untran_za = Normal(qzm_acc, qzv_acc.sqrt()).sample()
            z_acc = self.z_encoder_accessibility.z_transformation(untran_za)

            qzm_expr = qzm_expr.unsqueeze(0).expand(
                (n_samples, qzm_expr.size(0), qzm_expr.size(1))
            )
            qzv_expr = qzv_expr.unsqueeze(0).expand(
                (n_samples, qzv_expr.size(0), qzv_expr.size(1))
            )
            untran_zr = Normal(qzm_expr, qzv_expr.sqrt()).sample()
            z_expr = self.z_encoder_expression.z_transformation(untran_zr)

            libsize_expr = libsize_expr.unsqueeze(0).expand(
                (n_samples, libsize_expr.size(0), libsize_expr.size(1))
            )
            libsize_acc = libsize_acc.unsqueeze(0).expand(
                (n_samples, libsize_acc.size(0), libsize_acc.size(1))
            )

        ## Sample from the average distribution
        qzp_m = (qzm_acc + qzm_expr) / 2
        qzp_v = (qzv_acc + qzv_expr) / (2 ** 0.5)
        zp = Normal(qzp_m, qzp_v.sqrt()).rsample()

        ## choose the correct latent representation based on the modality
        qz_m = self._mix_modalities(qzp_m, qzm_expr, qzm_acc, mask_expr, mask_acc)
        qz_v = self._mix_modalities(qzp_v, qzv_expr, qzv_acc, mask_expr, mask_acc)
        z = self._mix_modalities(zp, z_expr, z_acc, mask_expr, mask_acc)

        outputs = dict(
            z=z,
            qz_m=qz_m,
            qz_v=qz_v,
            z_expr=z_expr,
            qzm_expr=qzm_expr,
            qzv_expr=qzv_expr,
            z_acc=z_acc,
            qzm_acc=qzm_acc,
            qzv_acc=qzv_acc,
            libsize_expr=libsize_expr,
            libsize_acc=libsize_acc,
        )
        return outputs

    def _get_generative_input(self, tensors, inference_outputs, transform_batch=None):
        z = inference_outputs["z"]
        qz_m = inference_outputs["qz_m"]
        libsize_expr = inference_outputs["libsize_expr"]
        labels = tensors[_CONSTANTS.LABELS_KEY]

        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        cont_key = _CONSTANTS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = _CONSTANTS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        input_dict = dict(
            z=z,
            qz_m=qz_m,
            batch_index=batch_index,
            cont_covs=cont_covs,
            cat_covs=cat_covs,
            libsize_expr=libsize_expr,
            labels=labels,
        )
        return input_dict

    @auto_move_data
    def generative(
        self,
        z,
        qz_m,
        batch_index,
        cont_covs=None,
        cat_covs=None,
        libsize_expr=None,
        labels=None,
        use_z_mean=False,
    ):
        """Runs the generative model."""
        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        latent = z if not use_z_mean else qz_m
        decoder_input = (
            latent if cont_covs is None else torch.cat([latent, cont_covs], dim=-1)
        )

        # Accessibility Decoder
        p = self.z_decoder_accessibility(decoder_input, batch_index, *categorical_input)

        # Expression Decoder
        px_scale, _, px_rate, px_dropout = self.z_decoder_expression(
            "gene", decoder_input, libsize_expr, batch_index, *categorical_input, labels
        )

        return dict(
            p=p,
            px_scale=px_scale,
            px_r=torch.exp(self.px_r),
            px_rate=px_rate,
            px_dropout=px_dropout,
        )

    def loss(
        self, tensors, inference_outputs, generative_outputs, kl_weight: float = 1.0
    ):
        # Get the data
        x = tensors[_CONSTANTS.X_KEY]

        x_rna = x[:, : self.n_input_genes]
        x_chr = x[:, self.n_input_genes :]

        mask_expr = x_rna.sum(dim=1) > 0
        mask_acc = x_chr.sum(dim=1) > 0

        # Compute Accessibility loss
        x_accessibility = x[:, self.n_input_genes :]
        p = generative_outputs["p"]
        libsize_acc = inference_outputs["libsize_acc"]
        rl_accessibility = self.get_reconstruction_loss_accessibility(
            x_accessibility, p, libsize_acc
        )

        # Compute Expression loss
        px_rate = generative_outputs["px_rate"]
        px_r = generative_outputs["px_r"]
        px_dropout = generative_outputs["px_dropout"]
        x_expression = x[:, : self.n_input_genes]
        rl_expression = self.get_reconstruction_loss_expression(
            x_expression, px_rate, px_r, px_dropout
        )

        # mix losses to get the correct loss for each cell
        recon_loss = self._mix_modalities(
            rl_accessibility + rl_expression,  # paired
            rl_expression,  # expression
            rl_accessibility,  # accessibility
            mask_expr,
            mask_acc,
        )

        # Compute KLD between Z and N(0,I)
        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]
        kl_div_z = kld(
            Normal(qz_m, torch.sqrt(qz_v)),
            Normal(0, 1),
        ).sum(dim=1)

        # Compute KLD between distributions for paired data
        qzm_expr = inference_outputs["qzm_expr"]
        qzv_expr = inference_outputs["qzv_expr"]
        qzm_acc = inference_outputs["qzm_acc"]
        qzv_acc = inference_outputs["qzv_acc"]
        kld_paired = kld(
            Normal(qzm_expr, torch.sqrt(qzv_expr)), Normal(qzm_acc, torch.sqrt(qzv_acc))
        ) + kld(
            Normal(qzm_acc, torch.sqrt(qzv_acc)), Normal(qzm_expr, torch.sqrt(qzv_expr))
        )
        kld_paired = torch.where(
            torch.logical_and(mask_acc, mask_expr),
            kld_paired.T,
            torch.zeros_like(kld_paired).T,
        ).sum(dim=0)

        # KL WARMUP
        kl_local_for_warmup = kl_div_z
        weighted_kl_local = kl_weight * kl_local_for_warmup

        # PENALTY
        # distance_penalty = kl_weight * torch.pow(z_acc - z_expr, 2).sum(dim=1)

        # TOTAL LOSS
        loss = torch.mean(recon_loss + weighted_kl_local + kld_paired)

        kl_local = dict(kl_divergence_z=kl_div_z)
        kl_global = torch.tensor(0.0)
        return LossRecorder(loss, recon_loss, kl_local, kl_global)

    def get_reconstruction_loss_expression(self, x, px_rate, px_r, px_dropout):
        rl = 0.0
        if self.gene_likelihood == "zinb":
            rl = (
                -ZeroInflatedNegativeBinomial(
                    mu=px_rate, theta=px_r, zi_logits=px_dropout
                )
                .log_prob(x)
                .sum(dim=-1)
            )
        elif self.gene_likelihood == "nb":
            rl = -NegativeBinomial(mu=px_rate, theta=px_r).log_prob(x).sum(dim=-1)
        elif self.gene_likelihood == "poisson":
            rl = -Poisson(px_rate).log_prob(x).sum(dim=-1)
        return rl

    def get_reconstruction_loss_accessibility(self, x, p, d):
        f = torch.sigmoid(self.region_factors) if self.region_factors is not None else 1
        return torch.nn.BCELoss(reduction="none")(p * d * f, (x > 0).float()).sum(
            dim=-1
        )

    @staticmethod
    def _mix_modalities(x_paired, x_expr, x_acc, mask_expr, mask_acc):
        """
        Mixes modality-specific vectors according to the modality masks.

        in positions where both `mask_expr` and `mask_acc` are True (corresponding to cell
        for which both expression and accessibility data is available), values from `x_paired`
        will be used. If only `mask_expr` is True, use values from `x_expr`, and if only
        `mask_acc` is True, use values from `x_acc`.

        Parameters
        ----------
        x_paired
            the values for paired cells (both modalities available), will be used in
            positions where both `mask_expr` and `mask_acc` are True.
        x_expr
            the values for expression-only cells, will be used in positions where
            only `mask_expr` is True.
        x_acc
            the values for accessibility-only cells, will be used on positions where
            only `mask_acc` is True.
        mask_expr
            the expression mask, indicating which cells have expression data
        mask_acc
            the accessibility mask, indicating which cells have accessibility data
        """
        x = torch.where(mask_expr.T, x_expr.T, x_acc.T).T
        x = torch.where(torch.logical_and(mask_acc, mask_expr), x_paired.T, x.T).T
        return x
