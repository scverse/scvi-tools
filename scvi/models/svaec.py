import torch
import torch.nn as nn
from torch.distributions import Normal, Multinomial, kl_divergence as kl

from scvi.metrics.log_likelihood import log_zinb_positive
from scvi.models.modules import Decoder, Encoder, Classifier, DecoderSCVI
from scvi.models.utils import enumerate_discrete, one_hot


class SVAEC(nn.Module):
    '''
    "Stacked" variational autoencoder for classification - SVAEC
    (from the stacked generative model M1 + M2)
    '''

    def __init__(self, n_input, n_labels, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, n_batch=0,
                 y_prior=None, use_cuda=False):
        super(SVAEC, self).__init__()
        self.n_labels = n_labels
        self.n_input = n_input

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
        self.register_buffer('px_r', torch.randn(n_input, ))

        # Classifier takes n_latent as input
        self.classifier = Classifier(n_latent, n_hidden, self.n_labels, n_layers, dropout_rate)
        self.encoder_z2_z1 = Encoder(n_input=n_latent, n_cat=self.n_labels, n_latent=n_latent, n_layers=n_layers)
        self.decoder_z1_z2 = Decoder(n_latent, n_latent, n_cat=self.n_labels, n_layers=n_layers)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.cuda()
            self.y_prior = self.y_prior.cuda()

    def classify(self, x):
        qz_m, _, z = self.z_encoder(x)
        return self.classifier(z)

    def sample_from_posterior_z(self, x, y=None):
        # Here we compute as little as possible to have q(z|x)
        qz_m, qz_v, z = self.z_encoder(x)
        return z

    def sample_from_posterior_l(self, x):
        # Here we compute as little as possible to have q(z|x)
        ql_m, ql_v, library = self.l_encoder.forward(x)
        return library

    def get_sample_scale(self, x, y=None, batch_index=None):
        z = self.sample_from_posterior_z(x, y)
        px = self.decoder.px_decoder(z, batch_index, y)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y=None, batch_index=None):
        z = self.sample_from_posterior_z(x)
        library = self.sample_from_posterior_l(x)
        px = self.decoder.px_decoder(z, batch_index, y)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        xs, ys = (x, y)
        xs_ = torch.log(1 + xs)
        qz1_m, qz1_v, z1_ = self.z_encoder(xs_)
        z1 = z1_
        # Enumerate choices of label
        if not is_labelled:
            ys = enumerate_discrete(xs, self.n_labels)
            xs = xs.repeat(self.n_labels, 1)
            if batch_index is not None:
                batch_index = batch_index.repeat(self.n_labels, 1)
            local_l_var = local_l_var.repeat(self.n_labels, 1)
            local_l_mean = local_l_mean.repeat(self.n_labels, 1)
            qz1_m = qz1_m.repeat(self.n_labels, 1)
            qz1_v = qz1_v.repeat(self.n_labels, 1)
            z1 = z1.repeat(self.n_labels, 1)
        else:
            ys = one_hot(ys, self.n_labels)

        xs_ = torch.log(1 + xs)

        qz2_m, qz2_v, z2 = self.encoder_z2_z1(z1, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)

        # Sampling
        ql_m, ql_v, library = self.l_encoder(xs_)  # let's keep that ind. of y

        px_scale, px_rate, px_dropout = self.decoder(self.dispersion, z1, library, batch_index)

        reconst_loss = -log_zinb_positive(xs, px_rate, torch.exp(self.px_r), px_dropout)

        # KL Divergence
        mean = torch.zeros_like(qz2_m)
        scale = torch.zeros_like(qz2_v)

        kl_divergence_z2 = kl(Normal(qz2_m, torch.sqrt(qz2_v)), Normal(mean, scale)).sum(dim=1)
        kl_divergence_z1 = kl(Normal(qz1_m, torch.sqrt(qz1_v)), Normal(pz1_m, torch.sqrt(pz1_v))).sum(dim=1)
        kl_divergence_l = kl(Normal(ql_m, torch.sqrt(ql_v)), Normal(local_l_mean, torch.sqrt(local_l_var))).sum(dim=1)
        kl_divergence = kl_divergence_z2 + kl_divergence_z1 + kl_divergence_l

        if is_labelled:
            return reconst_loss, kl_divergence

        reconst_loss = reconst_loss.view(self.n_labels, -1)
        kl_divergence = kl_divergence.view(self.n_labels, -1)

        probs = self.classifier(z1_)
        reconst_loss = (reconst_loss.t() * probs).sum(dim=1)
        kl_divergence = (kl_divergence.t() * probs).sum(dim=1)

        kl_divergence += kl(Multinomial(probs=probs), Multinomial(probs=self.y_prior))

        return reconst_loss, kl_divergence
