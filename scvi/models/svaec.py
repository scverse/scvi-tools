import torch
import torch.nn as nn
from torch.distributions import Normal, Multinomial, kl_divergence as kl

from scvi.metrics.log_likelihood import log_zinb_positive
from scvi.models.classifier import Classifier, LinearLogRegClassifier
from scvi.models.modules import Decoder, Encoder, DecoderSCVI
from scvi.models.utils import broadcast_labels
from .base import SemiSupervisedModel


class SVAEC(nn.Module, SemiSupervisedModel):
    '''
    "Stacked" variational autoencoder for classification - SVAEC
    (from the stacked generative model M1 + M2)
    '''

    def __init__(self, n_input, n_labels, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, n_batch=0,
                 y_prior=None, use_cuda=False, logreg_classifier=True):
        super(SVAEC, self).__init__()
        self.n_labels = n_labels
        self.n_input = n_input
        self.n_latent_layers = 2

        self.y_prior = y_prior if y_prior is not None else (1 / self.n_labels) * torch.ones(self.n_labels)
        # Automatically desactivate if useless
        self.n_batch = 0 if n_batch == 1 else n_batch
        self.z_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers,
                                 dropout_rate=dropout_rate)
        self.l_encoder = Encoder(n_input, n_hidden=n_hidden, n_latent=1, n_layers=1,
                                 dropout_rate=dropout_rate)
        self.decoder = DecoderSCVI(n_latent, n_input, n_hidden=n_hidden, n_layers=n_layers,
                                   dropout_rate=dropout_rate, n_batch=n_batch)

        self.dispersion = 'gene'
        self.px_r = torch.nn.Parameter(torch.randn(n_input, ))

        # Classifier takes n_latent as input
        if logreg_classifier:
            self.classifier = LinearLogRegClassifier(n_latent, self.n_labels)
        else:
            self.classifier = Classifier(n_latent, n_hidden, self.n_labels, n_layers, dropout_rate)

        self.encoder_z2_z1 = Encoder(n_input=n_latent, n_cat=self.n_labels, n_latent=n_latent, n_layers=n_layers)
        self.decoder_z1_z2 = Decoder(n_latent, n_latent, n_cat=self.n_labels, n_layers=n_layers)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.cuda()
            self.y_prior = self.y_prior.cuda()

    def classify(self, x):
        x_ = torch.log(1 + x)
        qz_m, _, z = self.z_encoder(x_)
        if not self.training:
            z = qz_m
        return self.classifier(z)

    def get_latents(self, x, y=None):
        x = torch.log(1 + x)
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder(x)
        if not self.training:
            z = qz_m
        return [z]

    def sample_from_posterior_z(self, x, y=None):
        x = torch.log(1 + x)
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder(x)
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
        z = self.sample_from_posterior_z(x)
        library = self.sample_from_posterior_l(x)
        px = self.decoder.px_decoder(z, batch_index, y)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        x_ = torch.log(1 + x)
        qz1_m, qz1_v, z1 = self.z_encoder(x_)
        ql_m, ql_v, library = self.l_encoder(x_)

        # Enumerate choices of label
        ys, z1s = (
            broadcast_labels(
                y, z1, n_broadcast=self.n_labels
            )
        )
        qz2_m, qz2_v, z2 = self.encoder_z2_z1(z1s, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)
        px_scale, px_rate, px_dropout = self.decoder(self.dispersion, z1, library, batch_index)
        reconst_loss = -log_zinb_positive(x, px_rate, torch.exp(self.px_r), px_dropout)

        # KL Divergence
        mean = torch.zeros_like(qz2_m)
        scale = torch.ones_like(qz2_v)

        kl_divergence_z2 = kl(Normal(qz2_m, torch.sqrt(qz2_v)), Normal(mean, scale)).sum(dim=1)
        loss_z1_unweight = - Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = Normal(qz1_m, torch.sqrt(qz1_v)).log_prob(z1).sum(dim=-1)
        kl_divergence_l = kl(Normal(ql_m, torch.sqrt(ql_v)), Normal(local_l_mean, torch.sqrt(local_l_var))).sum(dim=1)

        if is_labelled:
            return reconst_loss + loss_z1_weight + loss_z1_unweight, kl_divergence_z2 + kl_divergence_l

        probs = self.classifier(z1)
        reconst_loss += (loss_z1_weight + ((loss_z1_unweight).view(self.n_labels, -1).t() * probs).sum(dim=1))

        kl_divergence = (kl_divergence_z2.view(self.n_labels, -1).t() * probs).sum(dim=1)
        kl_divergence += kl(Multinomial(probs=probs), Multinomial(probs=self.y_prior))

        return reconst_loss, kl_divergence
