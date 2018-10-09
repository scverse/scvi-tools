# -*- coding: utf-8 -*-
"""Main module."""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, Poisson, LogNormal, kl_divergence as kl

from scvi.models.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.models.modules import Encoder, Decoder, DecoderSCVI
from scvi.models.utils import one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class VAECITE(nn.Module):
    r"""Variational auto-encoder model for CITE-seq data

    :param n_input: Number of input genes
    :param n_batch: Number of batches
    :param n_labels: Number of labels
    :param n_hidden: Number of nodes per hidden layer
    :param n_latent: Dimensionality of the latent space
    :param n_layers: Number of hidden layers used for encoder and decoder NNs
    :param dropout_rate: Dropout rate for neural networks
    :param rna_dispersion: One of the following

        * ``'gene'`` - rna_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - rna_dispersion can differ between different batches
        * ``'gene-label'`` - rna_dispersion can differ between different labels
        * ``'gene-cell'`` - rna_dispersion can differ for every gene in every cell

    :param log_variational: Log variational distribution
    :param reconstruction_loss:  One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels, use_cuda=True )

    """

    def __init__(self, n_input_genes, protein_indexes, n_batch=0, n_labels=0, n_hidden=128, n_latent=10,
                 n_layers=1, dropout_rate=0.1, rna_dispersion="gene", adt_dispersion='protein-cell',
                 log_variational=True, reconstruction_loss_umi="zinb",
                 reconstruction_loss_adt="poisson", model_library=True,
                 adt_mean_lib=None, adt_var_lib=None):
        super(VAECITE, self).__init__()
        self.rna_dispersion = rna_dispersion
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.reconstruction_loss_umi = reconstruction_loss_umi
        # Automatically deactivate if useless
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_input_genes = n_input_genes
        self.n_input_proteins = len(protein_indexes)
        self.reconstruction_loss_adt = reconstruction_loss_adt
        self.model_library = model_library
        self.adt_mean_lib = adt_mean_lib
        self.adt_var_lib = adt_var_lib
        self.adt_dispersion = adt_dispersion


        if self.rna_dispersion == "gene":
            self.px_r_umi = torch.nn.Parameter(torch.randn(n_input_genes, ))
        elif self.rna_dispersion == "gene-batch":
            self.px_r_umi = torch.nn.Parameter(torch.randn(n_input_genes, n_batch))
        elif self.rna_dispersion == "gene-label":
            self.px_r_umi = torch.nn.Parameter(torch.randn(n_input_genes, n_labels))
        else:  # gene-cell
            pass

        if self.adt_dispersion == "protein":
            self.px_r_adt = torch.nn.Parameter(torch.randn(self.n_input_proteins, ))
        elif self.adt_dispersion == "protein-batch":
            self.px_r_adt = torch.nn.Parameter(torch.randn(self.n_input_proteins, n_batch))
        elif self.adt_dispersion == "protein-label":
            self.px_r_adt = torch.nn.Parameter(torch.randn(self.n_input_proteins, n_labels))
        else:  # protein-cell
            # Decoder module expects this name
            self.adt_dispersion = 'gene-cell'

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        self.z_encoder = Encoder(n_input_genes + self.n_input_proteins, n_latent, n_layers=n_layers, n_hidden=n_hidden,
                                 dropout_rate=dropout_rate)
        self.l_umi_encoder = Encoder(n_input_genes, 1, n_hidden=n_hidden, n_layers=n_layers,
                                     dropout_rate=dropout_rate)
        self.l_adt_encoder = Encoder(self.n_input_proteins, 1, n_hidden=10, n_layers=n_layers,
                                     dropout_rate=dropout_rate)
        self.umi_decoder = DecoderSCVI(n_latent, n_input_genes, n_layers=n_layers, n_hidden=n_hidden,
                                       n_cat_list=[n_batch])
        if self.reconstruction_loss_adt == 'log_normal':
            self.adt_decoder = Decoder(n_latent, self.n_input_proteins, n_layers=n_layers, n_hidden=10, n_cat_list=[n_batch])
        else:
            self.adt_decoder = DecoderSCVI(n_latent, self.n_input_proteins, n_layers=n_layers, n_hidden=10, n_cat_list=[n_batch])

    def get_latents(self, x, y=None):
        r""" returns the result of ``sample_from_posterior_z`` inside a list

        :param x: tensor of values with shape ``(batch_size, n_input_genes + n_input_proteins)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: one element list of tensor
        :rtype: list of :py:class:`torch.Tensor`
        """
        return [self.sample_from_posterior_z(x, y)]

    def sample_from_posterior_z(self, x, y=None):
        r""" samples the tensor of latent values from the posterior
        #doesn't really sample, returns the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input_genes + n_input_proteins)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, n_latent)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
        qz_m, qz_v, z = self.z_encoder(x, y)  # y only used in VAEC
        return z

    def sample_from_posterior_l_umi(self, x):
        r""" samples the tensor of library sizes from the posterior
        #doesn't really sample, returns the tensor of the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, 1)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
        ql_m, ql_v, library = self.l_umi_encoder(x)
        return library

    def sample_from_posterior_l_adt(self, x):
        r""" samples the tensor of library sizes from the posterior
        #doesn't really sample, returns the tensor of the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input_adt)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, 1)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
        ql_m, ql_v, library = self.l_adt_encoder(x)
        return library

    def get_sample_scale_umi(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of predicted frequencies of expression

        :param x: tensor of values with shape ``(batch_size, n_input_genes)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of predicted frequencies of expression with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, y=y, n_samples=n_samples)[0]['umi']

    def get_sample_scale_adt(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of predicted frequencies of expression

        :param x: tensor of values with shape ``(batch_size, n_input_proteins)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of predicted frequencies of expression with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, y=y, n_samples=n_samples)[0]['adt']

    def get_sample_rate_umi(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of means of the negative binomial distribution

        :param x: tensor of values with shape ``(batch_size, n_input)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param n_samples: number of samples
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, y=y, n_samples=n_samples)[2]['umi']

    def get_sample_rate_adt(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of means of the negative binomial distribution

        :param x: tensor of values with shape ``(batch_size, n_input)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param n_samples: number of samples
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, y=y, n_samples=n_samples)[2]['adt']

    def _reconstruction_loss(self, x, px_rate, px_r, px_dropout, px_scale):
        # Reconstruction Loss
        umi = x[:, :self.n_input_genes]
        adt = x[:, self.n_input_genes:]
        if self.reconstruction_loss_umi == 'zinb':
            reconst_loss_umi = -log_zinb_positive(umi, px_rate['umi'], px_r['umi'], px_dropout['umi'])
        else:
            reconst_loss_umi = -log_nb_positive(umi, px_rate['umi'], px_r['umi'])

        if self.reconstruction_loss_adt == 'poisson':
            reconst_loss_adt = -torch.sum(Poisson(px_rate['adt']).log_prob(adt), dim=1)
        if self.reconstruction_loss_adt == 'nb':
            reconst_loss_adt = -log_nb_positive(adt, px_rate['adt'], px_r['adt'])
        if self.reconstruction_loss_adt == 'zinb':
            reconst_loss_adt = -log_zinb_positive(adt, px_rate['adt'], px_r['adt'], px_dropout['adt'])
        if self.reconstruction_loss_adt == 'log_normal':
            reconst_loss_adt = -torch.sum(LogNormal(px_rate['adt'], px_r['adt']).log_prob(adt), dim=1)

        return reconst_loss_umi, reconst_loss_adt

    def inference(self, x, batch_index=None, y=None, n_samples=1):
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        umi_ = x_[:, :self.n_input_genes]
        adt_ = x_[:, self.n_input_genes:]
        # Sampling - Encoder gets concatenated genes + proteins
        qz_m, qz_v, z = self.z_encoder(x_, y)
        ql_m = {}
        ql_v = {}
        ql_m['umi'], ql_v['umi'], library_umi = self.l_umi_encoder(umi_)
        if self.model_library:
            ql_m['adt'], ql_v['adt'], library_adt = self.l_adt_encoder(adt_)

        # TODO WHAT IS THIS CODE
        # if n_samples > 1:
        #     qz_m = qz_m.unsqueeze(0).expand(
        #         (n_samples, qz_m.size(0), qz_m.size(1)))
        #     qz_v = qz_v.unsqueeze(0).expand(
        #         (n_samples, qz_v.size(0), qz_v.size(1)))
        #     z = Normal(qz_m, qz_v.sqrt()).sample()
        #     ql_m_umi = ql_m_umi.unsqueeze(0).expand(
        #         (n_samples, ql_m_umi.size(0), ql_m_umi.size(1)))
        #     ql_v_umi = ql_v_umi.unsqueeze(0).expand(
        #         (n_samples, ql_v_umi.size(0), ql_v_umi.size(1)))
        #     library_umi = Normal(ql_m_umi, ql_v_umi.sqrt()).sample()
        #     if self.model_library:
        #         ql_m_adt = ql_m_adt.unsqueeze(0).expand(
        #             (n_samples, ql_m_adt.size(0), ql_m_adt.size(1)))
        #         ql_v_adt = ql_v_adt.unsqueeze(0).expand(
        #             (n_samples, ql_v_adt.size(0), ql_v_adt.size(1)))
        #         library_adt = Normal(ql_m_adt, ql_v_adt.sqrt()).sample()

        px_scale = {}
        px_r = {}
        px_rate = {}
        px_dropout = {}
        px_scale['umi'], px_r['umi'], px_rate['umi'], px_dropout['umi'] = self.umi_decoder(self.rna_dispersion, z, library_umi, batch_index, y)
        if self.rna_dispersion == "gene-label":
            px_r['umi'] = F.linear(one_hot(y, self.n_labels), self.px_r_umi)  # px_r gets transposed - last dimension is nb genes
        elif self.rna_dispersion == "gene-batch":
            px_r['umi'] = F.linear(one_hot(batch_index, self.n_batch), self.px_r_umi)
        elif self.rna_dispersion == "gene":
            px_r['umi'] = self.px_r_umi
        px_r['umi'] = torch.exp(px_r['umi'])

        if self.reconstruction_loss_adt != 'log_normal':
            px_scale['adt'], px_r['adt'], px_rate['adt'], px_dropout['adt'] = self.adt_decoder(self.adt_dispersion, z, library_adt, batch_index, y)
            if self.adt_dispersion == "protein-label":
                px_r['adt'] = F.linear(one_hot(y, self.n_labels), self.px_r_adt)  # px_r gets transposed - last dimension is nb genes
            elif self.adt_dispersion == "protein-batch":
                px_r['adt'] = F.linear(one_hot(batch_index, self.n_batch), self.px_r_adt)
            elif self.adt_dispersion == "protein":
                px_r['adt'] = self.px_r_adt
            px_r['adt'] = torch.exp(px_r['adt'])
        else:
            mean, var = self.adt_decoder(z, batch_index, y)
            px_rate['adt'] = mean + library_adt
            px_scale['adt'] = mean
            if self.adt_dispersion == "protein":
                px_r['adt'] = self.px_r_adt
                px_r['adt'] = torch.exp(px_r['adt'])
            else:
                px_r['adt'] = var

        return px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v

    def forward(self, x, local_l_mean_umi, local_l_var_umi, batch_index=None, y=None):
        r""" Returns the reconstruction loss and the Kullback divergences

        :param x: tensor of values with shape (batch_size, n_input)
        :param local_l_mean: tensor of means of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param local_l_var: tensor of variancess of the prior distribution of latent variable l
         with shape (batch_size, 1)
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape (batch_size, n_labels)
        :return: the reconstruction loss and the Kullback divergences
        :rtype: 2-tuple of :py:class:`torch.FloatTensor`
        """
        # Parameters for z latent distribution

        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v = self.inference(x, batch_index, y)
        reconst_loss_umi, reconst_loss_adt = self._reconstruction_loss(x, px_rate, px_r, px_dropout, px_scale)

        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)),
                             Normal(mean, scale)).sum(dim=1)
        kl_divergence_l_umi = kl(Normal(ql_m['umi'], torch.sqrt(ql_v['umi'])), Normal(local_l_mean_umi, torch.sqrt(local_l_var_umi))).sum(dim=1)
        if self.model_library:
            local_l_mean_adt = self.adt_mean_lib * torch.ones_like(ql_m['adt'])
            local_l_var_adt = self.adt_var_lib * torch.ones_like(ql_v['adt'])
            kl_divergence_l_adt = kl(Normal(ql_m['adt'], torch.sqrt(ql_v['adt'])), Normal(local_l_mean_adt, torch.sqrt(local_l_var_adt))).sum(dim=1)
        else:
            kl_divergence_l_adt = 0
        kl_divergence = kl_divergence_z

        return reconst_loss_umi + kl_divergence_l_umi, reconst_loss_adt + kl_divergence_l_adt, kl_divergence
