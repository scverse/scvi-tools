# -*- coding: utf-8 -*-
"""Main module."""

import torch
import torch.nn.functional as F
from torch.distributions import Normal, Poisson, kl_divergence as kl

from scvi.models import VAE
from scvi.models.classifier import Classifier
from scvi.models.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.models.modules import Encoder, DecoderSCVI
from scvi.models.utils import one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class VAEF(VAE):
    r"""Variational auto-encoder model.

    Args:
        :n_input: Number of input genes for scRNA-seq data.
        :n_input_fish: Number of input genes for smFISH data
        :n_batch: Default: ``0``.
        :n_labels: Default: ``0``.
        :n_hidden: Number of hidden. Default: ``128``.
        :n_latent: Default: ``1``.
        :n_layers: Number of layers. Default: ``1``.
        :dropout_rate: Default: ``0.1``.
        :dispersion: Default: ``"gene"``.
        :log_variational: Default: ``True``.
        :reconstruction_loss: Default: ``"zinb"``.
        :reconstruction_loss_fish: Default: ``"poisson"``.

    Examples:
        >>> gene_dataset_seq = CortexDataset()
        >>> gene_dataset_fish = SmfishDataset()
        >>> vae = VAE(gene_dataset_seq.nb_genes, gene_dataset_fish.nb_genes,
        ... n_labels=gene_dataset.n_labels, use_cuda=True )

    """

    def __init__(self, n_input, indexes_fish_train=None, n_batch=0, n_labels=0, n_hidden=128, n_latent=10,
                 n_layers=1, n_layers_decoder=1, dropout_rate=0.3,
                 dispersion="gene", log_variational=True, reconstruction_loss="zinb",
                 reconstruction_loss_fish="poisson", model_library=False):
        super().__init__(n_input, dispersion=dispersion, n_latent=n_hidden, n_hidden=n_hidden,
                         log_variational=log_variational, dropout_rate=dropout_rate, n_layers=1,
                         reconstruction_loss=reconstruction_loss, n_batch=n_batch, n_labels=n_labels)
        self.n_input = n_input
        self.n_input_fish = len(indexes_fish_train)
        self.indexes_to_keep = indexes_fish_train
        self.reconstruction_loss_fish = reconstruction_loss_fish
        self.model_library = model_library
        self.n_latent = n_latent
        # First layer of the encoder isn't shared
        self.z_encoder_fish = Encoder(self.n_input_fish, n_hidden, n_hidden=n_hidden, n_layers=1,
                                      dropout_rate=dropout_rate)
        # The last layers of the encoder are shared
        self.z_final_encoder = Encoder(n_hidden, n_latent, n_hidden=n_hidden, n_layers=n_layers,
                                       dropout_rate=dropout_rate)
        self.l_encoder_fish = Encoder(self.n_input_fish, 1, n_hidden=n_hidden, n_layers=1,
                                      dropout_rate=dropout_rate)
        self.l_encoder = Encoder(n_input, 1, n_hidden=n_hidden, n_layers=1,
                                 dropout_rate=dropout_rate)

        self.decoder = DecoderSCVI(n_latent, n_input, n_layers=n_layers_decoder, n_hidden=n_hidden,
                                   n_cat_list=[n_batch])

        self.classifier = Classifier(n_latent, n_labels=n_labels, n_hidden=128, n_layers=3)

    def get_latents(self, x, y=None):
        r""" returns the result of ``sample_from_posterior_z`` inside a list

        :param x: tensor of values with shape ``(batch_size, n_input)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: one element list of tensor
        :rtype: list of :py:class:`torch.Tensor`
        """
        return [self.sample_from_posterior_z(x, y)]

    def sample_from_posterior_z(self, x, y=None, mode="scRNA"):
        r""" samples the tensor of latent values from the posterior
        #doesn't really sample, returns the mean of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param mode: string that indicates the type of data we analyse
        :return: tensor of shape ``(batch_size, n_latent)``
        :rtype: :py:class:`torch.Tensor`
        """
        x = torch.log(1 + x)
        # First layer isn't shared
        if mode == "scRNA":
            z, _, _ = self.z_encoder(x)
        elif mode == "smFISH":
            z, _, _ = self.z_encoder_fish(x[:, self.indexes_to_keep])
        # The last layers of the encoder are shared
        qz_m, qz_v, z = self.z_final_encoder(z)
        if not self.training:
            z = qz_m
        return z

    def sample_from_posterior_l(self, x, mode="scRNA"):
        r""" samples the tensor of library sizes from the posterior
        #doesn't really sample, returns the tensor of the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param mode: string that indicates the type of data we analyse
        :return: tensor of shape ``(batch_size, 1)``
        :rtype: :py:class:`torch.Tensor`
        """
        x = torch.log(1 + x)
        if mode == "scRNA":
            ql_m, ql_v, library = self.l_encoder(x)
        elif mode == "smFISH":
            ql_m, ql_v, library = self.l_encoder_fish(x)
        return library

    def get_sample_scale(self, x, mode="scRNA", batch_index=None, y=None):
        r"""Returns the tensor of predicted frequencies of expression

        :param x: tensor of values with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param mode: string that indicates the type of data we analyse
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of predicted frequencies of expression with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :rtype: :py:class:`torch.Tensor`
        """
        z = self.sample_from_posterior_z(x, y, mode)  # y only used in VAEC
        px = self.decoder.px_decoder(z, batch_index, y)  # y only used in VAEC
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y=None, mode="scRNA"):
        r"""Returns the tensor of means of the negative binomial distribution

        :param x: tensor of values with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param mode: string that indicates the type of data we analyse
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :rtype: :py:class:`torch.Tensor`
        """
        if mode == "scRNA":
            library = torch.log(torch.sum(x, dim=1)).view(-1, 1)
            batch_index = torch.zeros_like(library)
        else:
            library = torch.log(torch.sum(x[:, self.indexes_to_keep], dim=1)).view(-1, 1)
            batch_index = torch.ones_like(library)

        if self.model_library:
            library = self.sample_from_posterior_l(x, mode=mode)
        px_scale = self.get_sample_scale(x, batch_index=batch_index, y=y)
        return px_scale * torch.exp(library)

    def get_sample_rate_fish(self, x, y=None):
        r"""Returns the tensor of means of the negative binomial distribution

        :param x: tensor of values with shape ``(batch_size, n_input_fish)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_input_fish)``
        :rtype: :py:class:`torch.Tensor`
        """
        library = torch.log(torch.sum(x[:, self.indexes_to_keep], dim=1)).view(-1, 1)
        batch_index = torch.ones_like(library)

        if self.model_library:
            library = self.sample_from_posterior_l(x, mode="smFISH")
        px_scale = self.get_sample_scale(x, mode="smFISH", batch_index=batch_index, y=y)
        px_scale = px_scale[:, self.indexes_to_keep] / torch.sum(px_scale[:, self.indexes_to_keep],
                                                                 dim=1).view(-1, 1)
        return px_scale * torch.exp(library)

    def classify(self, x, mode="scRNA"):
        r"""Classifies the cells based on their latent representation
        #for each cell, it gives the probability distribution over the different labels

        :param x: tensor of values with shape (batch_size, n_input)
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param mode: string that indicates the type of data we analyse
        :return: tensor of probabilities with shape``(batch_size, n_labels)``
        :rtype: :py:class:`torch.Tensor`
        """
        z = self.sample_from_posterior_z(x, mode)
        return self.classifier(z)

    def _reconstruction_loss(self, x, px_rate, px_r, px_dropout, batch_index, y, mode="scRNA", weighting=1):
        if self.dispersion == "gene-label":
            px_r = F.linear(one_hot(y, self.n_labels), self.px_r)  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r

        # Reconstruction Loss
        if mode == "scRNA":
            if self.reconstruction_loss == 'zinb':
                reconst_loss = -log_zinb_positive(x, px_rate, torch.exp(px_r), px_dropout)
            elif self.reconstruction_loss == 'nb':
                reconst_loss = - log_nb_positive(x, px_rate, torch.exp(px_r))

        else:
            if self.reconstruction_loss_fish == 'poisson':
                reconst_loss = -torch.sum(Poisson(px_rate).log_prob(x), dim=1)
            elif self.reconstruction_loss_fish == 'gaussian':
                reconst_loss = -torch.sum(Normal(px_rate, 10).log_prob(x), dim=1)
        return reconst_loss

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None, mode="scRNA", weighting=1):
        r""" Returns the reconstruction loss and the Kullback divergences

        :param x: tensor of values with shape ``(batch_size, n_input)``
        or ``(batch_size, n_input_fish)`` depending on the mode
        :param local_l_mean: tensor of means of the prior distribution of latent variable l
        with shape (batch_size, 1)
        :param local_l_var: tensor of variances of the prior distribution of latent variable l
        with shape (batch_size, 1)
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape (batch_size, n_labels)
        :param mode: string that indicates the type of data we analyse
        :param weighting: used in none of these methods
        :return: the reconstruction loss and the Kullback divergences
        :rtype: 2-tuple of :py:class:`torch.FloatTensor`
        """
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)
        # Sampling
        if mode == "scRNA":
            qz_m, qz_v, z = self.z_encoder(x_)
            library = torch.log(torch.sum(x, dim=1)).view(-1, 1)
            batch_index = torch.zeros_like(library)
        if mode == "smFISH":
            qz_m, qz_v, z = self.z_encoder_fish(x_[:, self.indexes_to_keep])
            library = torch.log(torch.sum(x[:, self.indexes_to_keep], dim=1)).view(-1, 1)
            batch_index = torch.ones_like(library)
        if self.model_library:
            if mode == "scRNA":
                ql_m, ql_v, library = self.l_encoder(x_)
            elif mode == "smFISH":
                ql_m, ql_v, library = self.l_encoder_fish(x_[:, self.indexes_to_keep])

        qz_m, qz_v, z = self.z_final_encoder(z)
        px_scale, px_r, px_rate, px_dropout = self.decoder(self.dispersion, z, library, batch_index)

        # rescaling the expected frequencies
        if mode == "smFISH":
            if self.model_library:
                px_rate = px_scale[:, self.indexes_to_keep] * torch.exp(library)
                reconst_loss = self._reconstruction_loss(x[:, self.indexes_to_keep], px_rate, px_r, px_dropout,
                                                         batch_index, y, mode)
            else:
                px_scale = px_scale[:, self.indexes_to_keep] / torch.sum(
                    px_scale[:, self.indexes_to_keep], dim=1).view(-1, 1)
                px_rate = px_scale * torch.exp(library)
                reconst_loss = self._reconstruction_loss(x[:, self.indexes_to_keep], px_rate, px_r, px_dropout,
                                                         batch_index, y, mode)

        else:
            reconst_loss = self._reconstruction_loss(x, px_rate, px_r, px_dropout, batch_index, y, mode, weighting)

        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(dim=1)
        if self.model_library:
            kl_divergence_l = kl(Normal(ql_m, torch.sqrt(ql_v)), Normal(local_l_mean,
                                                                        torch.sqrt(local_l_var))).sum(dim=1)
            kl_divergence = kl_divergence_z + kl_divergence_l
        else:
            kl_divergence = kl_divergence_z

        return reconst_loss, kl_divergence
