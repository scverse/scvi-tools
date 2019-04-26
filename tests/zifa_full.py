import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, kl_divergence as kl

from scvi.models.log_likelihood import log_zinb_positive, log_nb_positive
from scvi.models.modules import Encoder, DecoderSCVI
from scvi.models.utils import one_hot
import numpy as np
import matplotlib.pyplot as plt

import collections
from typing import Iterable, Tuple

import torch
from torch import nn as nn
from torch.distributions import Normal

class FCLayers(nn.Module):
    r"""A helper class to build fully-connected layers for a neural network.

    :param n_in: The dimensionality of the input
    :param n_out: The dimensionality of the output
    :param n_layers: The number of fully-connected hidden layers
    :param n_hidden: The number of nodes per hidden layer
    :param dropout_rate: Dropout rate to apply to each of the hidden layers
    """

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        use_batch_norm=True,
    ):
        super().__init__()
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]

        self.fc_layers = nn.Sequential(
            collections.OrderedDict(
                [
                    (
                        "Layer {}".format(i),
                        nn.Sequential(
                            nn.Linear(n_in, n_out),
                            nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001)
                            if use_batch_norm
                            else None,
                            nn.ReLU(),
                            nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None,
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(
                        zip(layers_dim[:-1], layers_dim[1:])
                    )
                ]
            )
        )

    def forward(self, x: torch.Tensor):
        r"""Forward computation on ``x``.
        
        :param x: tensor of values with shape ``(n_in,)``
        :return: tensor of shape ``(n_out,)``
        :rtype: :py:class:`torch.Tensor`
        """

        for layers in self.fc_layers:
            for layer in layers:
                if layer is not None:
                    if isinstance(layer, nn.BatchNorm1d):
                        if x.dim() == 3:
                            x = torch.cat(
                                [(layer(slice_x)).unsqueeze(0) for slice_x in x], dim=0
                            )
                        else:
                            x = layer(x)
                    else:
                        x = layer(x)
        return x


# Encoder
class Encoder(nn.Module):
    r"""Encodes data of ``n_input`` dimensions into a latent space of ``n_output``
    dimensions using a fully-connected neural network of ``n_hidden`` layers.

    :param n_input: The dimensionality of the input (data space)
    :param n_output: The dimensionality of the output (latent space)
    :param n_layers: The number of fully-connected hidden layers
    :param n_hidden: The number of nodes per hidden layer
    :dropout_rate: Dropout rate to apply to each of the hidden layers
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
    ):
        super().__init__()

        self.encoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
        )
        self.mean_encoder = nn.Linear(n_hidden, n_output)
        self.var_encoder = nn.Linear(n_hidden, n_output)

    def reparameterize(self, mu, var):
        return Normal(mu, var.sqrt()).rsample()

    def forward(self, x: torch.Tensor, *cat_list: int):
        r"""The forward computation for a single sample.

         #. Encodes the data into latent space using the encoder network
         #. Generates a mean \\( q_m \\) and variance \\( q_v \\) (clamped to \\( [-5, 5] \\))
         #. Samples a new value from an i.i.d. multivariate normal \\( \\sim N(q_m, \\mathbf{I}q_v) \\)

        :param x: tensor with shape (n_input,)
        :return: tensors of shape ``(n_latent,)`` for mean and var, and sample
        :rtype: 3-tuple of :py:class:`torch.Tensor`
        """

        # Parameters for latent distribution
        q = self.encoder(x)
        q_m = self.mean_encoder(q)
        q_v = torch.exp(
            self.var_encoder(q)
        )  # (computational stability safeguard)torch.clamp(, -5, 5)
        latent = self.reparameterize(q_m, q_v)
        return q_m, q_v, latent


class DecoderSCVI(nn.Module):
    r"""Decodes data from latent space of ``n_input`` dimensions ``n_output``
    dimensions using a fully-connected neural network of ``n_hidden`` layers.

    :param n_input: The dimensionality of the input (latent space)
    :param n_output: The dimensionality of the output (data space)
    :param n_layers: The number of fully-connected hidden layers
    :param n_hidden: The number of nodes per hidden layer
    :param dropout_rate: Dropout rate to apply to each of the hidden layers
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_layers: int = 1,
        n_hidden: int = 128,
        decay_mode: str = 'gene',
    ):
        super().__init__()
        self.px_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
        )

        # mean gamma
        self.px_scale_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output), nn.Softmax(dim=-1)
        )

        # dispersion: here we only deal with gene-cell dispersion case
        self.px_r_decoder = nn.Linear(n_hidden, n_output)

        # dropout
        # self.px_dropout_decoder = nn.Linear(n_hidden, n_output)
        
        self.alpha = nn.Sequential(
            nn.Linear(n_hidden, n_output), nn.Sigmoid()
        )
        
        if decay_mode == 'gene':
            self.decay = torch.nn.Parameter(torch.randn(n_output,))
            self.decay_bias = torch.nn.Parameter(torch.randn(n_output,))
        else:
            self.decay = torch.nn.Parameter(torch.randn(1))
            self.decay_bias = torch.nn.Parameter(torch.randn(1))

    def forward(
        self, dispersion: str, z: torch.Tensor, library: torch.Tensor, 
    ):
        r"""The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns parameters for the ZINB distribution of expression
         #. If ``dispersion != 'gene-cell'`` then value for that param will be ``None``

        :param dispersion: One of the following

            * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
            * ``'gene-batch'`` - dispersion can differ between different batches
            * ``'gene-label'`` - dispersion can differ between different labels
            * ``'gene-cell'`` - dispersion can differ for every gene in every cell

        :param z: tensor with shape ``(n_input,)``
        :param library: library size
        :return: parameters for the ZINB distribution of expression
        :rtype: 4-tuple of :py:class:`torch.Tensor`
        """

        # The decoder returns values for the parameters of the ZINB distribution
        px = self.px_decoder(z)
        px_scale = self.px_scale_decoder(px)
        # Clamp to high value: exp(12) ~ 160000 to avoid nans (computational stability)
        px_rate = torch.exp(library) * px_scale  # torch.clamp( , max=12)
        
        # px_dropout = self.px_dropout_decoder(px)
        px_dropout_unlogit = torch.exp(- torch.exp(self.decay_bias) - torch.exp(self.decay) * torch.pow(px_rate, 2))
        px_dropout = torch.log( px_dropout_unlogit / (1-px_dropout_unlogit + 1e-8) + 1e-8)
        
        
        
        px_r = self.px_r_decoder(px) if dispersion == "gene-cell" else None
        # alpha = self.alpha(px)
        return px_scale, px_r, px_rate, px_dropout#, alpha

class VAE(nn.Module):
    r"""Variational auto-encoder model.

    :param n_genes: Number of input genes
    :param n_batch: Number of batches
    :param n_labels: Number of labels
    :param n_hidden: Number of nodes per hidden layer
    :param n_latent: Dimensionality of the latent space
    :param n_layers: Number of hidden layers used for encoder and decoder NNs
    :param dropout_rate: Dropout rate for neural networks
    :param dispersion: One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell

    :param log_variational: Log variational distribution
    
    """

    def __init__(self, n_genes: int, n_batch: int = 0, n_labels: int = 0,
                 n_hidden: int = 128, n_latent: int = 10, n_layers: int = 1,
                 dropout_rate: float = 0.1, dispersion: str = "gene",
                 log_variational: bool = True, reconstruction_loss: str = "zinb", decay_mode:str = 'gene'):
        super().__init__()
        self.dispersion = dispersion
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.reconstruction_loss = reconstruction_loss
        # Automatically deactivate if useless
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_latent_layers = 1  # not sure what this is for, no usages?

        if self.dispersion == "gene":
            self.px_r = torch.nn.Parameter(torch.randn(n_genes, ))
        elif self.dispersion == "gene-batch":
            self.px_r = torch.nn.Parameter(torch.randn(n_genes, n_batch))
        elif self.dispersion == "gene-label":
            self.px_r = torch.nn.Parameter(torch.randn(n_genes, n_labels))
        else:  # gene-cell
            pass
        
        #self.alpha_p = torch.nn.Parameter(torch.randn(n_genes, ))
        #self.alpha_s = nn.Sigmoid()(self.alpha_p)
        
        # z encoder goes from the n_genes-dimensional data to an n_latent-d
        # latent space representation
        self.z_encoder = Encoder(n_genes, n_latent, n_layers=n_layers, n_hidden=n_hidden,
                                 dropout_rate=dropout_rate)
        # l encoder goes from n_genes-dimensional data to 1-d library size
        self.l_encoder = Encoder(n_genes, 1, n_layers=1, n_hidden=n_hidden, dropout_rate=dropout_rate)
    
        # decoder goes from n_latent-dimensional space to n_genes-d data
        self.decoder = DecoderSCVI(n_latent, n_genes, n_layers=n_layers, n_hidden=n_hidden, decay_mode=decay_mode)
        
        # self.location_layer = FCLayers(2, n_genes, n_hidden=16, n_layers=3, dropout_rate=0.1)
        # alpha
        
        

    def get_latents(self, x, y=None):
        r""" returns the result of ``sample_from_posterior_z`` inside a list

        :param x: tensor of values with shape ``(batch_size, n_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: one element list of tensor
        :rtype: list of :py:class:`torch.Tensor`
        """
        return [self.sample_from_posterior_z(x)]

    def sample_from_posterior_z(self, x, y=None, give_mean=False):
        r""" samples the tensor of latent values from the posterior
        #doesn't really sample, returns the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param give_mean: is True when we want the mean of the posterior  distribution rather than sampling
        :return: tensor of shape ``(batch_size, n_latent)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
        qz_m, qz_v, z = self.z_encoder(x)  # y only used in VAEC
        if give_mean:
            z = qz_m
        return z

    def sample_from_posterior_l(self, x):
        r""" samples the tensor of library sizes from the posterior
        #doesn't really sample, returns the tensor of the means of the posterior distribution

        :param x: tensor of values with shape ``(batch_size, n_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :return: tensor of shape ``(batch_size, 1)``
        :rtype: :py:class:`torch.Tensor`
        """
        if self.log_variational:
            x = torch.log(1 + x)
        ql_m, ql_v, library = self.l_encoder(x)
        return library

    def get_sample_scale(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of predicted frequencies of expression

        :param x: tensor of values with shape ``(batch_size, n_genes)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param n_samples: number of samples
        :return: tensor of predicted frequencies of expression with shape ``(batch_size, n_genes)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, n_samples=n_samples)[0]

    def get_sample_rate(self, x, batch_index=None, y=None, n_samples=1):
        r"""Returns the tensor of means of the negative binomial distribution

        :param x: tensor of values with shape ``(batch_size, n_genes)``
        :param y: tensor of cell-types labels with shape ``(batch_size, n_labels)``
        :param batch_index: array that indicates which batch the cells belong to with shape ``batch_size``
        :param n_samples: number of samples
        :return: tensor of means of the negative binomial distribution with shape ``(batch_size, n_genes)``
        :rtype: :py:class:`torch.Tensor`
        """
        return self.inference(x, batch_index=batch_index, n_samples=n_samples)[2]

    def _reconstruction_loss(self, x, px_rate, px_r, px_dropout):
        # Reconstruction Loss
        if self.reconstruction_loss == 'zinb':
            reconst_loss = -log_zinb_positive(x, px_rate, px_r, px_dropout)
        elif self.reconstruction_loss == 'nb':
            reconst_loss = -log_nb_positive(x, px_rate, px_r)
        return reconst_loss

    def scale_from_z(self, sample_batch, fixed_batch):
        if self.log_variational:
            sample_batch = torch.log(1 + sample_batch)
        qz_m, qz_v, z = self.z_encoder(sample_batch)
        batch_index = fixed_batch * torch.ones_like(sample_batch[:, [0]])
        library = 4. * torch.ones_like(sample_batch[:, [0]])
        px_scale, _, _, _ = self.decoder('gene', z, library, batch_index)
        return px_scale

    def inference(self, x, batch_index=None, y=None, n_samples=1):
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        # Sampling
        qz_m, qz_v, z = self.z_encoder(x_)
        ql_m, ql_v, library = self.l_encoder(x_)

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            z = Normal(qz_m, qz_v.sqrt()).sample()
            ql_m = ql_m.unsqueeze(0).expand((n_samples, ql_m.size(0), ql_m.size(1)))
            ql_v = ql_v.unsqueeze(0).expand((n_samples, ql_v.size(0), ql_v.size(1)))
            library = Normal(ql_m, ql_v.sqrt()).sample()

        # px_scale, px_r, px_rate, px_dropout, alpha = self.decoder(self.dispersion, z, library)
        px_scale, px_r, px_rate, px_dropout = self.decoder(self.dispersion, z, library)
        if self.dispersion == "gene-label":
            px_r = F.linear(one_hot(y, self.n_labels), self.px_r)  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r
        px_r = torch.exp(px_r)

        return px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library #, alpha

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        r""" Returns the reconstruction loss and the Kullback divergences

        :param x: tensor of values with shape (batch_size, n_genes)
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

        # px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library, alpha = self.inference(x, batch_index)
        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library = self.inference(x, batch_index)
        
        # location_r = self.location_layer(locations)
        # px_scale = self.alpha_s*px_scale + (1-self.alpha_s)*location_r
    
        # px_r = px_r * location_r
        # px_scale = alpha*px_scale + (1-alpha)*location_r
        
        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(dim=1)
        kl_divergence_l = kl(Normal(ql_m, torch.sqrt(ql_v)), Normal(local_l_mean, torch.sqrt(local_l_var))).sum(dim=1)
        kl_divergence = kl_divergence_z

        reconst_loss = self._reconstruction_loss(x, px_rate, px_r, px_dropout)

        return reconst_loss + kl_divergence_l, kl_divergence

    
from scvi.inference import Trainer
class UnsupervisedTrainer(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\*\*kwargs: Other keywords arguments from the general Trainer class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, test_size=None, kl=None, **kwargs):
        super().__init__(model, gene_dataset, **kwargs)
        self.kl = kl
        if type(self) is UnsupervisedTrainer:
            self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size, test_size)
            self.train_set.to_monitor = ['ll']
            self.test_set.to_monitor = ['ll']

    @property
    def posteriors_loop(self):
        return ['train_set']

    def loss(self, tensors):
        # sample_batch, local_l_mean, local_l_var, batch_index, _, locations = tensors
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        # reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, locations, batch_index)
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / 400)  # self.n_epochs)
        
        
        

