from typing import Dict, Iterable, Optional, Union

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from torch.distributions import Normal, Poisson
from torch.distributions import kl_divergence as kld

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.distributions import (
    NegativeBinomial,
    ZeroInflatedNegativeBinomial,
    NegativeBinomialMixture,
)
from scvi.module._peakvae import Decoder as DecoderPeakVI
from scvi.module.base import BaseModuleClass, LossRecorder, auto_move_data
from scvi.nn import DecoderSCVI, Encoder, FCLayers, one_hot


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


class DecoderADT(torch.nn.Module):
    """
    Decoder for just surface proteins (ADT)
    """

    def __init__(
        self,
        n_input: int,
        n_output_proteins: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 2,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
        deep_inject_covariates: bool = False,
    ):
        super().__init__()
        self.n_output_proteins = n_output_proteins

        linear_args = dict(
            n_layers=1,
            use_activation=False,
            use_batch_norm=False,
            use_layer_norm=False,
            dropout_rate=0,
        )

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

        self.py_background_decoder = FCLayers(
            n_in=n_hidden + n_input,
            n_out=n_output_proteins,
            n_cat_list=n_cat_list,
            **linear_args,
        )

        # dropout (mixture component for proteins, ZI probability for genes)
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

        # background mean parameters second decoder
        self.py_back_mean_log_alpha = FCLayers(
            n_in=n_hidden + n_input,
            n_out=n_output_proteins,
            n_cat_list=n_cat_list,
            **linear_args,
        )
        self.py_back_mean_log_beta = FCLayers(
            n_in=n_hidden + n_input,
            n_out=n_output_proteins,
            n_cat_list=n_cat_list,
            **linear_args,
        )

        # background mean first decoder
        self.py_back_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
        )

    def forward(self, z: torch.Tensor, *cat_list: int):
        # z is the latent repr
        py_ = {}

        py_back = self.py_back_decoder(z, *cat_list)
        py_back_cat_z = torch.cat([py_back, z], dim=-1)

        py_["back_alpha"] = self.py_back_mean_log_alpha(py_back_cat_z, *cat_list)
        py_["back_beta"] = torch.exp(
            self.py_back_mean_log_beta(py_back_cat_z, *cat_list)
        )
        log_pro_back_mean = Normal(py_["back_alpha"], py_["back_beta"]).rsample()
        py_["rate_back"] = torch.exp(log_pro_back_mean)

        py_fore = self.py_fore_decoder(z, *cat_list)
        py_fore_cat_z = torch.cat([py_fore, z], dim=-1)
        py_["fore_scale"] = (
            self.py_fore_scale_decoder(py_fore_cat_z, *cat_list) + 1 + 1e-8
        )
        py_["rate_fore"] = py_["rate_back"] * py_["fore_scale"]

        p_mixing = self.sigmoid_decoder(z, *cat_list)
        p_mixing_cat_z = torch.cat([p_mixing, z], dim=-1)
        py_["mixing"] = self.py_background_decoder(p_mixing_cat_z, *cat_list)

        protein_mixing = 1 / (1 + torch.exp(-py_["mixing"]))
        py_["scale"] = torch.nn.functional.normalize(
            (1 - protein_mixing) * py_["rate_fore"], p=1, dim=-1
        )

        return py_, log_pro_back_mean


class MULTIVAE(BaseModuleClass):
    """
    Variational auto-encoder model for joint paired + unpaired RNA-seq, ATAC-seq, and CITE-seq data.

    Parameters
    ----------
    n_input_regions
        Number of input regions.
    n_input_genes
        Number of input genes.
    n_input_proteins
        Number of input proteins.
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

    def __init__(
        self,
        n_input_regions: int = 0,
        n_input_genes: int = 0,
        n_input_proteins: int = 0,
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
        protein_batch_mask: Dict[Union[str, int], np.ndarray] = None,
        protein_background_prior_mean: Optional[np.ndarray] = None,
        protein_background_prior_scale: Optional[np.ndarray] = None,
        protein_dispersion: str = "protein",
    ):
        super().__init__()

        # ####################################################################
        # INIT PARAMS
        self.n_input_regions = n_input_regions
        self.n_input_genes = n_input_genes
        self.n_input_proteins = n_input_proteins
        self.n_hidden = (
            int(np.sqrt(self.n_input_regions + self.n_input_genes))
            if n_hidden is None
            else n_hidden
        )
        self.n_batch = n_batch
        self.n_labels = n_labels

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
        encoder_cat_list = cat_list if encode_covariates else None
        n_input_decoder = self.n_latent + self.n_continuous_cov
        # ####################################################################

        # ####################################################################
        # RNA
        # expression encoder
        if self.n_input_genes == 0:
            input_exp = 1
        else:
            input_exp = self.n_input_genes
        n_input_encoder_exp = input_exp + n_continuous_cov * encode_covariates

        self.gene_likelihood = gene_likelihood
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
            n_input_decoder,
            n_input_genes,
            n_cat_list=cat_list,
            n_layers=n_layers_decoder,
            n_hidden=self.n_hidden,
            inject_covariates=self.deeply_inject_covariates,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
        )

        # expression dispersion parameters
        self.px_r = torch.nn.Parameter(torch.randn(n_input_genes))

        # expression library size encoder
        self.l_encoder_expression = LibrarySizeEncoder(
            n_input_encoder_exp,
            n_cat_list=encoder_cat_list,
            n_layers=self.n_layers_encoder,
            n_hidden=self.n_hidden,
            use_batch_norm=self.use_batch_norm_encoder,
            use_layer_norm=self.use_layer_norm_encoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )
        # ####################################################################

        # ####################################################################
        # ACCESSIBILITY
        # accessibility encoder
        if self.n_input_regions == 0:
            input_acc = 1
        else:
            input_acc = self.n_input_regions
        n_input_encoder_acc = input_acc + n_continuous_cov * encode_covariates

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

        # accessibility decoder
        self.z_decoder_accessibility = DecoderPeakVI(
            n_input=n_input_decoder,
            n_output=n_input_regions,
            n_hidden=self.n_hidden,
            n_cat_list=cat_list,
            n_layers=self.n_layers_decoder,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )

        # accessibility region-specific factors
        self.region_factors = None
        if region_factors:
            self.region_factors = torch.nn.Parameter(torch.zeros(input_acc))

        # accessibility library size encoder
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
        # ####################################################################

        # ####################################################################
        # PROTEIN
        self.protein_dispersion = protein_dispersion
        self.protein_batch_mask = protein_batch_mask
        if protein_background_prior_mean is None:
            if n_batch > 0:
                self.background_pro_alpha = torch.nn.Parameter(
                    torch.randn(n_input_proteins, n_batch)
                )
                self.background_pro_log_beta = torch.nn.Parameter(
                    torch.clamp(torch.randn(n_input_proteins, n_batch), -10, 1)
                )
            else:
                self.background_pro_alpha = torch.nn.Parameter(
                    torch.randn(n_input_proteins)
                )
                self.background_pro_log_beta = torch.nn.Parameter(
                    torch.clamp(torch.randn(n_input_proteins), -10, 1)
                )
        else:
            if protein_background_prior_mean.shape[1] == 1 and n_batch != 1:
                init_mean = protein_background_prior_mean.ravel()
                init_scale = protein_background_prior_scale.ravel()
            else:
                init_mean = protein_background_prior_mean
                init_scale = protein_background_prior_scale
            self.background_pro_alpha = torch.nn.Parameter(
                torch.from_numpy(init_mean.astype(np.float32))
            )
            self.background_pro_log_beta = torch.nn.Parameter(
                torch.log(torch.from_numpy(init_scale.astype(np.float32)))
            )

        ## protein encoder
        if self.n_input_proteins == 0:
            input_pro = 1
        else:
            input_pro = self.n_input_proteins
        n_input_encoder_pro = input_pro + n_continuous_cov * encode_covariates
        self.z_encoder_protein = Encoder(
            n_input=n_input_encoder_pro,
            n_layers=self.n_layers_encoder
            if n_input_encoder_pro < 1000
            else self.n_layers_encoder + 1,
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

        # protein decoder
        self.z_decoder_pro = DecoderADT(
            n_input=n_input_decoder,
            n_output_proteins=n_input_proteins,
            n_hidden=self.n_hidden,
            n_cat_list=cat_list,
            n_layers=self.n_layers_decoder,
            use_batch_norm=self.use_batch_norm_decoder,
            use_layer_norm=self.use_layer_norm_decoder,
            deep_inject_covariates=self.deeply_inject_covariates,
        )

        # protein dispersion parameters
        if self.protein_dispersion == "protein":
            self.py_r = torch.nn.Parameter(2 * torch.rand(self.n_input_proteins))
        elif self.protein_dispersion == "protein-batch":
            self.py_r = torch.nn.Parameter(2 * torch.rand(input_pro, n_batch))
        elif self.protein_dispersion == "protein-label":
            self.py_r = torch.nn.Parameter(2 * torch.rand(input_pro, n_labels))
        else:  # protein-cell
            pass
        # ####################################################################

    def _get_inference_input(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        if self.n_input_proteins == 0:
            y = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        cont_covs = tensors.get(_CONSTANTS.CONT_COVS_KEY)
        cat_covs = tensors.get(_CONSTANTS.CAT_COVS_KEY)
        input_dict = dict(
            x=x,
            y=y,
            batch_index=batch_index,
            cont_covs=cont_covs,
            cat_covs=cat_covs,
        )
        return input_dict

    @auto_move_data
    def inference(
        self,
        x,
        y,
        batch_index,
        cont_covs,
        cat_covs,
        n_samples=1,
    ) -> Dict[str, torch.Tensor]:

        # Get Data and Additional Covs
        if self.n_input_genes == 0:
            x_rna = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            x_rna = x[:, : self.n_input_genes]
        if self.n_input_regions == 0:
            x_chr = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            x_chr = x[:, self.n_input_genes :]

        mask_expr = x_rna.sum(dim=1) > 0
        mask_acc = x_chr.sum(dim=1) > 0
        mask_pro = y.sum(dim=1) > 0

        if cont_covs is not None and self.encode_covariates:
            encoder_input_expression = torch.cat((x_rna, cont_covs), dim=-1)
            encoder_input_accessibility = torch.cat((x_chr, cont_covs), dim=-1)
            encoder_input_protein = torch.cat((y, cont_covs), dim=-1)
        else:
            encoder_input_expression = x_rna
            encoder_input_accessibility = x_chr
            encoder_input_protein = y

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
        qzm_pro, qzv_pro, z_pro = self.z_encoder_protein(
            encoder_input_protein, batch_index, *categorical_input
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

            qzm_pro = qzm_pro.unsqueeze(0).expand(
                (n_samples, qzm_pro.size(0), qzm_pro.size(1))
            )
            qzv_pro = qzv_pro.unsqueeze(0).expand(
                (n_samples, qzv_pro.size(0), qzv_pro.size(1))
            )
            untran_zp = Normal(qzm_pro, qzv_pro.sqrt()).sample()
            z_pro = self.z_encoder_protein.z_transformation(untran_zp)

            libsize_expr = libsize_expr.unsqueeze(0).expand(
                (n_samples, libsize_expr.size(0), libsize_expr.size(1))
            )
            libsize_acc = libsize_acc.unsqueeze(0).expand(
                (n_samples, libsize_acc.size(0), libsize_acc.size(1))
            )

        ## Sample from the average distribution
        qzp_m123 = (qzm_acc + qzm_expr + qzm_pro) / 3
        qzp_v123 = (qzv_acc + qzv_expr + qzv_pro) / (3 ** 0.5)
        zp123 = Normal(qzp_m123, qzp_v123.sqrt()).rsample()

        qzp_m12 = (qzm_acc + qzm_expr) / 2
        qzp_v12 = (qzv_acc + qzv_expr) / (2 ** 0.5)
        zp12 = Normal(qzp_m12, qzp_v12.sqrt()).rsample()

        qzp_m13 = (qzm_acc + qzm_pro) / 2
        qzp_v13 = (qzv_acc + qzv_pro) / (2 ** 0.5)
        zp13 = Normal(qzp_m13, qzp_v13.sqrt()).rsample()

        qzp_m23 = (qzm_expr + qzm_pro) / 2
        qzp_v23 = (qzv_expr + qzv_pro) / (2 ** 0.5)
        zp23 = Normal(qzp_m23, qzp_v23.sqrt()).rsample()

        ## choose the correct latent representation based on the modality
        qz_m = self._mix_modalities123(
            qzp_m123,
            qzp_m12,
            qzp_m13,
            qzp_m23,
            qzm_expr,
            qzm_acc,
            qzm_pro,
            mask_expr,
            mask_acc,
            mask_pro,
        )

        qz_v = self._mix_modalities123(
            qzp_v123,
            qzp_v12,
            qzp_v13,
            qzp_v23,
            qzv_expr,
            qzv_acc,
            qzv_pro,
            mask_expr,
            mask_acc,
            mask_pro,
        )
        z = self._mix_modalities123(
            zp123, zp12, zp13, zp23, z_expr, z_acc, z_pro, mask_expr, mask_acc, mask_pro
        )

        if self.n_batch > 0:
            py_back_alpha_prior = F.linear(
                one_hot(batch_index, self.n_batch), self.background_pro_alpha
            )
            py_back_beta_prior = F.linear(
                one_hot(batch_index, self.n_batch),
                torch.exp(self.background_pro_log_beta),
            )
        else:
            py_back_alpha_prior = self.background_pro_alpha
            py_back_beta_prior = torch.exp(self.background_pro_log_beta)
        self.back_mean_prior = Normal(py_back_alpha_prior, py_back_beta_prior)

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
            z_pro=z_pro,
            qzm_pro=qzm_pro,
            qzv_pro=qzv_pro,
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

        # Protein Decoder
        py_, log_pro_back_mean = self.z_decoder_pro(
            decoder_input, batch_index, *categorical_input
        )

        py_["r"] = torch.exp(self.py_r)

        return dict(
            p=p,
            px_scale=px_scale,
            px_r=torch.exp(self.px_r),
            px_rate=px_rate,
            px_dropout=px_dropout,
            py_=py_,
            log_pro_back_mean=log_pro_back_mean,
        )

    def loss(
        self, tensors, inference_outputs, generative_outputs, kl_weight: float = 1.0
    ):
        # Get the data
        x = tensors[_CONSTANTS.X_KEY]
        if self.n_input_genes == 0:
            x_rna = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            x_rna = x[:, : self.n_input_genes]
        if self.n_input_regions == 0:
            x_chr = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            x_chr = x[:, self.n_input_genes :]
        if self.n_input_proteins == 0:
            y = torch.zeros(x.shape[0], 1, device=x.device, requires_grad=False)
        else:
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]

        mask_expr = x_rna.sum(dim=1) > 0
        mask_acc = x_chr.sum(dim=1) > 0
        mask_pro = y.sum(dim=1) > 0

        pro_batch_mask_minibatch = None
        if self.protein_batch_mask is not None:
            pro_batch_mask_minibatch = torch.zeros_like(y)
            for b in torch.unique(batch_index):
                b_indices = (batch_index == b).reshape(-1)
                pro_batch_mask_minibatch[b_indices] = torch.tensor(
                    self.protein_batch_mask[b.item()].astype(np.float32),
                    device=y.device,
                )

        # Compute Accessibility loss
        if self.n_input_regions == 0:
            rl_accessibility = torch.zeros(x.shape[0])
        else:
            x_accessibility = x[:, self.n_input_genes :]
            p = generative_outputs["p"]
            libsize_acc = inference_outputs["libsize_acc"]
            rl_accessibility = self.get_reconstruction_loss_accessibility(
                x_accessibility, p, libsize_acc
            )

        # Compute Expression loss
        if self.n_input_genes == 0:
            rl_expression = torch.zeros(x.shape[0])
        else:
            px_rate = generative_outputs["px_rate"]
            px_r = generative_outputs["px_r"]
            px_dropout = generative_outputs["px_dropout"]
            x_expression = x[:, : self.n_input_genes]
            rl_expression = self.get_reconstruction_loss_expression(
                x_expression, px_rate, px_r, px_dropout
            )

        # Compute Protein loss
        if self.n_input_proteins == 0:
            rl_protein = torch.zeros(x.shape[0])
        else:
            py_ = generative_outputs["py_"]
            rl_protein = self.get_reconstruction_loss_protein(
                y, py_, pro_batch_mask_minibatch
            )

        # mix losses to get the correct loss for each cell
        recon_loss = self._mix_modalities123(
            rl_accessibility + rl_expression + rl_protein,
            rl_accessibility + rl_expression,
            rl_accessibility + rl_protein,
            rl_expression + rl_protein,
            rl_expression,
            rl_accessibility,
            rl_protein,
            mask_expr,
            mask_acc,
            mask_pro,
        )

        # Compute KLD between Z and N(0,I)
        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]
        kl_div_z = kld(Normal(qz_m, torch.sqrt(qz_v)), Normal(0, 1)).sum(dim=1)

        # Compute KLD between distributions for paired data
        qzm_expr = inference_outputs["qzm_expr"]
        qzv_expr = inference_outputs["qzv_expr"]
        qzm_acc = inference_outputs["qzm_acc"]
        qzv_acc = inference_outputs["qzv_acc"]
        qzm_pro = inference_outputs["qzm_pro"]
        qzv_pro = inference_outputs["qzv_pro"]

        def symKLd(qzm1, qzv1, qzm2, qzv2):
            out = kld(
                Normal(qzm1, torch.sqrt(qzv1)), Normal(qzm2, torch.sqrt(qzv2))
            ) + kld(Normal(qzm2, torch.sqrt(qzv2)), Normal(qzm1, torch.sqrt(qzv1)))

            return out

        """
        kld_paired = (
            symKLd(qzm_expr, qzv_expr, qzm_acc, qzv_acc)
            + symKLd(qzm_expr, qzv_expr, qzm_pro, qzv_pro)
            + symKLd(qzm_acc, qzv_acc, qzm_pro, qzv_pro)
        )

        mask = torch.logical_or(
            torch.logical_or(
                torch.logical_or(
                    torch.logical_and(mask_acc, mask_expr),
                    torch.logical_and(mask_expr, mask_pro),
                ),
                torch.logical_and(mask_acc, mask_pro),
            ),
            torch.logical_and(torch.logical_and(mask_expr, mask_acc), mask_pro),
        )

        kld_paired = torch.where(
            mask,
            kld_paired.T,
            torch.zeros_like(kld_paired).T,
        ).sum(dim=0)
        """
        symKLd12 = symKLd(qzm_expr, qzv_expr, qzm_acc, qzv_acc)
        symKLd13 = symKLd(qzm_expr, qzv_expr, qzm_pro, qzv_pro)
        symKLd23 = symKLd(qzm_acc, qzv_acc, qzm_pro, qzv_pro)

        kld_paired = self._mix_modalities123(
            symKLd12 + symKLd13 + symKLd23,
            symKLd12,
            symKLd23,
            symKLd13,
            torch.zeros(x.shape[0], device=x.device, requires_grad=False),
            torch.zeros(x.shape[0], device=x.device, requires_grad=False),
            torch.zeros(x.shape[0], device=x.device, requires_grad=False),
            mask_expr,
            mask_acc,
            mask_pro,
        ).sum(dim=1)

        # KL WARMUP
        kl_local_for_warmup = kl_div_z
        weighted_kl_local = kl_weight * kl_local_for_warmup

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

    def get_reconstruction_loss_protein(self, y, py_, pro_batch_mask_minibatch):

        py_conditional = NegativeBinomialMixture(
            mu1=py_["rate_back"],
            mu2=py_["rate_fore"],
            theta1=py_["r"],
            mixture_logits=py_["mixing"],
        )

        reconst_loss_protein_full = -py_conditional.log_prob(y)
        if pro_batch_mask_minibatch is not None:
            temp_pro_loss_full = torch.zeros_like(reconst_loss_protein_full)
            temp_pro_loss_full.masked_scatter_(
                pro_batch_mask_minibatch.bool(), reconst_loss_protein_full
            )
            rl_protein = temp_pro_loss_full.sum(dim=-1)
        else:
            rl_protein = reconst_loss_protein_full.sum(dim=-1)

        return rl_protein

    @staticmethod
    def _mix_modalities123(
        x_123, x_12, x_13, x_23, x_1, x_2, x_3, mask_1, mask_2, mask_3
    ):
        """
        Mixes modality-specific vectors according to the modality masks.

        In positions where both `mask_1` and `mask_2` are True (corresponding to cell
        for which both modality data is available), values from `x_12` will be used, etc.
        If only `mask_1` is True, use values from `x_1`, and if only `mask_2` is True, use values from `x_2`.

        Parameters
        ----------
        x_123
            the values for fully paired cells (all modalities available), will be used in
            positions where both `mask_expr` and `mask_acc` are True.
        x_12
            the values for partially paired cells, will be used in
            positions where both `mask_expr` and `mask_acc` are True.
        x_13
            the values for partially paired cells, will be used in
            positions where both `mask_expr` and `mask_prot` are True.
        x_23
            the values for paired cells (both modalities available), will be used in
            positions where both `mask_prot` and `mask_acc` are True.
        x_1
            the values for modality1-only cells, will be used in positions where
            only `mask_1` is True.
        x_2
            the values for modality2-only cells, will be used on positions where
            only `mask_2` is True.
        x_3
            the values for modality3-only cells, will be used on positions where
            only `mask_3` is True.
        mask_1
            the mask for modality 1, indicating which cells have modality 1 data
        mask_2
            the mask for modality 2, indicating which cells have modality 2 data
        mask_3
            the mask for modality 3, indicating which cells have modality 3 data
        """

        # Single Modality
        x = torch.where(mask_1.T, x_1.T, x_2.T).T
        x = torch.where(mask_3.T, x_3.T, x.T).T

        # Double Modality
        mask_12 = torch.logical_and(mask_2, mask_1)
        x = torch.where(mask_12.T, x_12.T, x.T).T
        mask_13 = torch.logical_and(mask_3, mask_1)
        x = torch.where(mask_13.T, x_13.T, x.T).T
        mask_23 = torch.logical_and(mask_3, mask_2)
        x = torch.where(mask_23.T, x_23.T, x.T).T

        # Fully Paired
        mask_123 = torch.logical_and(mask_12, mask_3)
        x = torch.where(mask_123, x_123.T, x.T).T

        return x
