import torch
import torch.nn as nn
from torch.distributions import Normal, Multinomial, kl_divergence as kl

from scvi.metrics.log_likelihood import log_zinb_positive
from scvi.models.classifier import Classifier
from scvi.models.modules import Encoder, DecoderSCVI
from scvi.models.utils import broadcast_labels
from .base import SemiSupervisedModel


# VAE model - for classification: VAEC
class VAEC(nn.Module, SemiSupervisedModel):
    def __init__(self, n_input, n_labels, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, dispersion="gene",
                 log_variational=True, reconstruction_loss="zinb", n_batch=0, y_prior=None, use_cuda=False):
        super(VAEC, self).__init__()
        self.dispersion = dispersion
        self.log_variational = log_variational
        self.reconstruction_loss = reconstruction_loss
        self.n_latent_layers = 1
        # Automatically desactivate if useless
        self.n_batch = 0 if n_batch == 1 else n_batch
        self.n_labels = 0 if n_labels == 1 else n_labels
        if self.n_labels == 0:
            raise ValueError("VAEC is only implemented for > 1 label dataset")

        if self.dispersion == "gene":
            self.px_r = torch.nn.Parameter(torch.randn(n_input, ))

        self.z_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                                 dropout_rate=dropout_rate, n_cat=n_labels)
        self.l_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=1, n_layers=1,
                                 dropout_rate=dropout_rate)
        self.decoder = DecoderSCVI(n_latent, n_input, n_hidden=n_hidden, n_layers=n_layers,
                                   dropout_rate=dropout_rate, n_batch=n_batch, n_labels=n_labels)

        self.y_prior = y_prior if y_prior is not None else (1 / n_labels) * torch.ones(n_labels)
        self.classifier = Classifier(n_input, n_hidden, n_labels, n_layers=n_layers, dropout_rate=dropout_rate)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.cuda()
            self.y_prior = self.y_prior.cuda()

    def classify(self, x):
        x = torch.log(1 + x)
        return self.classifier(x)

    def get_latents(self, x, y=None):
        x = torch.log(1 + x)
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder(x, y)
        if self.training:
            z = qz_m
        return [z]

    def sample_from_posterior_z(self, x, y):
        x = torch.log(1 + x)
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder(x, y)
        return z

    def sample_from_posterior_l(self, x):
        x = torch.log(1 + x)
        # Here we compute as little as possible to have q(z|x)
        ql_m, ql_v, library = self.l_encoder(x)
        return library

    def get_sample_scale(self, x, y=None, batch_index=None):
        x = torch.log(1 + x)
        z = self.sample_from_posterior_z(x, y)
        px = self.decoder.px_decoder(z, batch_index, y)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y=None, batch_index=None):
        x = torch.log(1 + x)
        z = self.sample_from_posterior_z(x, y)
        library = self.sample_from_posterior_l(x)
        px = self.decoder.px_decoder(z, batch_index, y)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        # Prepare for sampling
        x_ = torch.log(1 + x)
        ql_m, ql_v, library = self.l_encoder(x_)

        # Enumerate choices of label
        ys, xs, library_s, batch_index_s = (
            broadcast_labels(
                y, x, library, batch_index, n_broadcast=self.n_labels
            )
        )

        if self.log_variational:
            xs_ = torch.log(1 + xs)

        # Sampling
        qz_m, qz_v, zs = self.z_encoder(xs_, ys)

        px_scale, px_rate, px_dropout = self.decoder(self.dispersion, zs, library_s, batch_index_s, y=ys)

        reconst_loss = -log_zinb_positive(xs, px_rate, torch.exp(self.px_r), px_dropout)

        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(dim=1)
        kl_divergence_l = kl(Normal(ql_m, torch.sqrt(ql_v)), Normal(local_l_mean, torch.sqrt(local_l_var))).sum(dim=1)

        if is_labelled:
            return reconst_loss, kl_divergence_z + kl_divergence_l

        reconst_loss = reconst_loss.view(self.n_labels, -1)

        probs = self.classifier(x_)
        reconst_loss = (reconst_loss.t() * probs).sum(dim=1)

        kl_divergence = (kl_divergence_z.view(self.n_labels, -1).t() * probs).sum(dim=1)
        kl_divergence += kl(Multinomial(probs=probs), Multinomial(probs=self.y_prior))
        kl_divergence += kl_divergence_l

        return reconst_loss, kl_divergence
