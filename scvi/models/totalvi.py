# -*- coding: utf-8 -*-
"""Main module."""
from typing import Dict, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, Poisson, kl_divergence as kl

from scvi.models.log_likelihood import (
    log_zinb_positive,
    log_nb_positive,
    log_mixture_nb,
)
from scvi.models.modules import DecoderTOTALVI, EncoderTOTALVI
from scvi.models.utils import one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class TOTALVI(nn.Module):
    r"""Latent variable model for CITE-seq data using auto-encoding Variational Bayes

    :param n_input_genes: Number of input genes
    :param protein_indexes: List of indexes (columns) which correspond to protein
    :param n_batch: Number of batches
    :param n_labels: Number of labels
    :param n_hidden: Number of nodes per hidden layer for the z encoder (protein+genes),
                     genes library encoder, z->genes+proteins decoder
    :param n_latent: Dimensionality of the latent space
    :param n_layers: Number of hidden layers used for encoder and decoder NNs
    :param dropout_rate: Dropout rate for neural networks
    :param genes_dispersion: One of the following

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels

    :param protein_dispersion: One of the following

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED

    :param log_variational: Log(data+1) prior to encoding for numerical stability. Not normalization.
    :param reconstruction_loss_genes:  One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution

    :param reconstruction_loss_protein:  One of

        * ``'mixture_nb'`` - Negative binomial mixture distribution

    :param latent_distribution:  One of

        * ``'normal'`` - Isotropic normal
        * ``'gsm'`` - Gaussian softmax N(0, 1) passed through softmax

    Examples:
        >>> dataset = Dataset10X(dataset_name="pbmc_10k_protein_v3", save_path=save_path)
        >>> totalvae = totalVI(gene_dataset.nb_genes, len(dataset.protein_names), use_cuda=True )
    """

    def __init__(
        self,
        n_input_genes: int,
        n_input_proteins: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 256,
        n_latent: int = 20,
        n_layers: int = 1,
        dropout_rate_decoder: float = 0.2,
        dropout_rate_encoder: float = 0.2,
        gene_dispersion: str = "gene",
        protein_dispersion: str = "protein",
        log_variational: bool = True,
        reconstruction_loss_gene: str = "nb",
        reconstruction_loss_protein: str = "mixture_nb",
        latent_distribution: str = "gsm",
    ):
        super().__init__()
        self.gene_dispersion = gene_dispersion
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.reconstruction_loss_gene = reconstruction_loss_gene
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_input_genes = n_input_genes
        self.n_input_proteins = n_input_proteins
        self.reconstruction_loss_protein = reconstruction_loss_protein
        self.protein_dispersion = protein_dispersion
        self.latent_distribution = latent_distribution

        if n_batch > 0:
            self.background_log_alpha = torch.nn.Parameter(
                torch.randn(n_input_proteins, n_batch)
            )
            self.background_log_beta = torch.nn.Parameter(
                torch.randn(n_input_proteins, n_batch)
            )
        else:
            self.background_log_alpha = torch.nn.Parameter(
                torch.randn(n_input_proteins)
            )
            self.background_log_beta = torch.nn.Parameter(torch.randn(n_input_proteins))

        if self.gene_dispersion == "gene":
            self.px_r_gene = torch.nn.Parameter(torch.randn(n_input_genes))
        elif self.gene_dispersion == "gene-batch":
            self.px_r_gene = torch.nn.Parameter(torch.randn(n_input_genes, n_batch))
        elif self.gene_dispersion == "gene-label":
            self.px_r_gene = torch.nn.Parameter(torch.randn(n_input_genes, n_labels))
        else:  # gene-cell
            pass

        if self.protein_dispersion == "protein":
            self.px_r_protein = torch.nn.Parameter(torch.ones(self.n_input_proteins))
        elif self.protein_dispersion == "protein-batch":
            self.px_r_protein = torch.nn.Parameter(
                torch.ones(self.n_input_proteins, n_batch)
            )
        elif self.protein_dispersion == "protein-label":
            self.px_r_protein = torch.nn.Parameter(
                torch.ones(self.n_input_proteins, n_labels)
            )
        else:  # protein-cell
            pass

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        self.encoder = EncoderTOTALVI(
            n_input_genes + self.n_input_proteins,
            n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_encoder,
            distribution=latent_distribution,
        )
        self.decoder = DecoderTOTALVI(
            n_latent,
            n_input_genes,
            self.n_input_proteins,
            n_layers=1,
            n_cat_list=[n_batch, n_labels],
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_decoder,
        )

    def get_latents(self, x: torch.Tensor, batch_index: Optional[torch.Tensor] = None):
        r""" returns the result of ``sample_from_posterior_z`` inside a list

        :param x: tensor of values with shape ``(batch_size, n_input_genes + n_input_proteins)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: one element list of tensor
        :rtype: list of :py:class:`torch.Tensor`
        """
        return [self.sample_from_posterior_z(x, batch_index)]

    def sample_from_posterior_z(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        give_mean: bool = False,
        n_samples: int = 5000,
    ) -> torch.Tensor:
        """ samples the tensor of latent values from the posterior
        #doesn't really sample, returns the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input_genes + n_input_proteins)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, n_latent)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
            y = torch.log(1 + y)
        qz_m, qz_v, z, _, _, _, _ = self.encoder(torch.cat((x, y), dim=-1), batch_index)
        if give_mean:
            if self.latent_distribution == "gsm":
                samples = Normal(qz_m, qz_v.sqrt()).sample([n_samples])
                z = self.encoder.transformation(samples)
                z = z.mean(dim=0)
            else:
                z = qz_m
        return z

    def sample_from_posterior_l(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        give_mean: bool = True,
    ) -> torch.Tensor:
        r""" samples the tensor of library size from the posterior
        #doesn't really sample, returns the tensor of the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, 1)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
            y = torch.log(1 + y)
        _, _, _, _, ql_m, ql_v, library = self.encoder(
            torch.cat((x, y), dim=-1), batch_index
        )
        if give_mean is True:
            return torch.exp(ql_m + 0.5 * ql_v)
        else:
            return library

    def get_sample_scale(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> torch.Tensor:
        r"""Returns the tensor of predicted frequencies of expression for RNA/Proteins

        Protein data is concatenated horizontally to gene data

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of values with shape ``(batch_size, n_input_proteins)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param label: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of predicted frequencies of expression with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        px_scale = self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )["px_scale"]
        px_scale = torch.cat((px_scale["gene"], px_scale["protein"]), dim=-1)
        return px_scale

    def get_background_mean(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> torch.Tensor:
        px_rate = self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )["px_rate"]
        px_rate = px_rate["background"]
        return px_rate

    def get_sample_dropout(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> torch.Tensor:
        r"""Returns the tensor of mixing components for RNA/proteins

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of values with shape ``(batch_size, n_input_proteins)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param label: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of mixing probabilities on logits scale ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`

        The first n_input_genes columns of the return value are reserved for ZINB mixing components,
        while the last n_input_proteins columns are reserved for the background probability of NB mixture.
        """

        px_dropout = self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )["px_dropout"]
        px_dropout = torch.cat((px_dropout["gene"], px_dropout["protein"]), dim=-1)
        return px_dropout

    def get_sample_rate(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> torch.Tensor:
        """Returns the tensor of means of the negative binomial distribution

        Protein means are horizontally concatenated to gene means.

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of values with shape ``(batch_size, n_input_proteins)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param label: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        px_rate = self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )["px_rate"]
        px_rate = torch.cat((px_rate["gene"], px_rate["protein"]), dim=-1)
        return px_rate

    def get_sample_dispersion(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
        mode: str = "protein",
    ) -> torch.Tensor:
        """Returns the tensor of dispersions/variances depending on the model

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of values with shape ``(batch_size, n_input_proteins)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param label: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :param mode: one of "protein" or "gene"
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )["px_r"][mode]

    def scale_from_z(
        self, x: torch.Tensor, y: torch.Tensor, fixed_batch: torch.Tensor
    ) -> torch.Tensor:
        if self.log_variational:
            x = torch.log(1 + x)
            y = torch.log(1 + y)
        qz_m, qz_v, z, _, _, _, _ = self.encoder(torch.cat((x, y), dim=-1))
        batch_index = fixed_batch * torch.ones_like(x[:, [0]])
        library = 4.0 * torch.ones_like(x[:, [0]])
        px_scale = self.decoder(z, library, batch_index)[0]
        px_scale = torch.cat((px_scale["gene"], px_scale["protein"]), dim=-1)
        return px_scale

    def get_reconstruction_loss(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        px_rate: Dict[str, torch.Tensor],
        px_r: Dict[str, torch.Tensor],
        px_dropout: Dict[str, torch.Tensor],
        px_scale: Dict[str, torch.Tensor],
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        # Reconstruction Loss
        gene = x
        protein = y
        if self.reconstruction_loss_gene == "zinb":
            reconst_loss_gene = -log_zinb_positive(
                gene, px_rate["gene"], px_r["gene"], px_dropout["gene"]
            ).sum(dim=-1)
        else:
            reconst_loss_gene = -log_nb_positive(
                gene, px_rate["gene"], px_r["gene"]
            ).sum(dim=-1)

        if self.reconstruction_loss_protein == "poisson":
            reconst_loss_protein = (
                -Poisson(px_rate["protein"]).log_prob(protein).sum(dim=-1)
            )
        if self.reconstruction_loss_protein == "nb":
            reconst_loss_protein = -log_nb_positive(
                protein, px_rate["protein"], px_r["protein"]
            ).sum(dim=-1)
        if self.reconstruction_loss_protein == "mixture_nb":

            reconst_loss_protein = -log_mixture_nb(
                protein,
                px_rate["background"],
                px_rate["protein"],
                px_r["protein"],
                px_r["protein"],
                px_dropout["protein"],
            ).sum(dim=-1)

        return reconst_loss_gene, reconst_loss_protein

    def inference(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples=1,
    ) -> Dict[str, Union[torch.Tensor, Dict[str, torch.Tensor]]]:
        """ Internal helper function to compute necessary inference quanitities

        We use the variable name `px_...` to refer to parameters of the respective distributions
        (gene, protein). The rate refers to the mean of the NB, dropout refers to Bernoulli
        mixing parameters. `scale` refers to the quanity upon which differential expression is
        performed. For genes, this can be viewed as the mean of the underlying gamma distribution.
        The `px_...` variables are dictionaries, typically with `background, protein, gene` keys.
        """
        x_ = x
        y_ = y
        if self.log_variational:
            x_ = torch.log(1 + x_)
            y_ = torch.log(1 + y_)

        px_r = {}

        # Sampling - Encoder gets concatenated genes + proteins
        qz_m, qz_v, z, untran_z, ql_m, ql_v, library_gene = self.encoder(
            torch.cat((x_, y_), dim=-1), batch_index
        )

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            untran_z = Normal(qz_m, qz_v.sqrt()).sample()
            if self.latent_distribution != "normal":
                z = self.encoder.transformation(untran_z)
            else:
                z = untran_z
            ql_m = ql_m.unsqueeze(0).expand((n_samples, ql_m.size(0), ql_m.size(1)))
            ql_v = ql_v.unsqueeze(0).expand((n_samples, ql_v.size(0), ql_v.size(1)))
            library_gene = Normal(ql_m, ql_v.sqrt()).sample()

        if self.gene_dispersion == "gene-label":
            # px_r gets transposed - last dimension is nb genes
            px_r["gene"] = F.linear(one_hot(label, self.n_labels), self.px_r_gene)
        elif self.gene_dispersion == "gene-batch":
            px_r["gene"] = F.linear(one_hot(batch_index, self.n_batch), self.px_r_gene)
        elif self.gene_dispersion == "gene":
            px_r["gene"] = self.px_r_gene
        px_r["gene"] = torch.exp(px_r["gene"])

        if self.protein_dispersion == "protein-label":
            # px_r gets transposed - last dimension is nb genes
            px_r["protein"] = F.linear(one_hot(label, self.n_labels), self.px_r_protein)
        elif self.protein_dispersion == "protein-batch":
            px_r["protein"] = F.linear(
                one_hot(batch_index, self.n_batch), self.px_r_protein
            )
        elif self.protein_dispersion == "protein":
            px_r["protein"] = self.px_r_protein

        px_r["protein"] = torch.exp(px_r["protein"])

        # Background regularization
        if self.n_batch > 0:
            back_alpha = F.linear(
                one_hot(batch_index, self.n_batch), self.background_log_alpha
            )
            back_beta = F.linear(
                one_hot(batch_index, self.n_batch), torch.exp(self.background_log_beta)
            )
        else:
            back_alpha = self.background_log_alpha
            back_beta = torch.exp(self.background_log_beta)
        self.back_mean_prior = Normal(back_alpha, back_beta)

        px_scale, px_rate, px_dropout, py_alpha, py_beta, log_back_mean = self.decoder(
            z, library_gene, batch_index, label
        )

        protein_mixing = 1 / (1 + torch.exp(-px_dropout["protein"]))
        px_scale["protein"] = F.normalize(
            (1 - protein_mixing) * px_rate["protein"], p=1, dim=-1
        )

        return dict(
            px_scale=px_scale,
            px_r=px_r,
            px_rate=px_rate,
            px_dropout=px_dropout,
            qz_m=qz_m,
            qz_v=qz_v,
            z=z,
            untran_z=untran_z,
            ql_m=ql_m,
            ql_v=ql_v,
            library_gene=library_gene,
            py_alpha=py_alpha,
            py_beta=py_beta,
            log_back_mean=log_back_mean,
        )

    def forward(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        local_l_mean_gene: torch.Tensor,
        local_l_var_gene: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
    ):
        r""" Returns the reconstruction loss and the Kullback divergences

        :param x: tensor of values with shape (batch_size, n_input_genes)
        :param y: tensor of values with shape (batch_size, n_input_proteins)
        :param local_l_mean_gene: tensor of means of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param local_l_var_gene: tensor of variancess of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param label: tensor of cell-types labels with shape (batch_size, n_labels)
        :return: the reconstruction loss and the Kullback divergences
        :rtype: 4-tuple of :py:class:`torch.FloatTensor`
        """
        # Parameters for z latent distribution

        outputs = self.inference(x, y, batch_index, label)
        qz_m = outputs["qz_m"]
        qz_v = outputs["qz_v"]
        ql_m = outputs["ql_m"]
        ql_v = outputs["ql_v"]
        px_rate = outputs["px_rate"]
        px_scale = outputs["px_scale"]
        px_r = outputs["px_r"]
        px_dropout = outputs["px_dropout"]

        reconst_loss_gene, reconst_loss_protein = self.get_reconstruction_loss(
            x, y, px_rate, px_r, px_dropout, px_scale
        )

        # KL Divergence
        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(0, 1)).sum(dim=1)
        kl_divergence_l_gene = kl(
            Normal(ql_m, torch.sqrt(ql_v)),
            Normal(local_l_mean_gene, torch.sqrt(local_l_var_gene)),
        ).sum(dim=1)

        back_reg = kl(
            Normal(outputs["py_alpha"], outputs["py_beta"]), self.back_mean_prior
        ).sum(dim=-1)

        kl_local = kl_divergence_z + kl_divergence_l_gene
        return (reconst_loss_gene, reconst_loss_protein, kl_local, back_reg)
