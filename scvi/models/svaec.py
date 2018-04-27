import torch

from scvi.log_likelihood import log_zinb_positive
from scvi.modules import Decoder, Encoder, Classifier
from scvi.utils import enumerate_discrete, one_hot
from . import VAE


class SVAEC(VAE):
    '''
    "Stacked" variational autoencoder for classification - SVAEC
    (from the stacked generative model M1 + M2)
    '''

    def __init__(self, n_input, n_labels=7, y_prior=None, n_latent=10, n_layers=3, **kwargs):
        super(SVAEC, self).__init__(n_input, n_latent=n_latent, n_layers=n_layers, **kwargs)
        self.n_labels = n_labels
        self.n_input = n_input

        self.y_prior = y_prior if y_prior is not None else (1 / self.n_labels) * torch.ones(self.n_labels)

        # Classifier takes n_latent as input
        self.classifier = Classifier(n_latent, self.n_hidden, self.n_labels, self.n_layers, self.dropout_rate)
        self.encoder_z2_z1 = Encoder(n_input=n_latent + self.n_labels, n_latent=n_latent, n_layers=n_layers)
        self.decoder_z1_z2 = Decoder(n_latent + self.n_labels, n_latent, n_layers=n_layers)

    def classify(self, x):
        qz_m, _, z = self.z_encoder(x)
        return self.classifier(z)

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

        if self.batch and batch_index is None:
            raise ("This VAE was trained to take batches into account:"
                   "please provide batch indexes when running the forward pass")

        qz2_m, qz2_v, z2 = self.encoder_z2_z1(torch.cat((z1, ys), 1))
        pz1_m, pz1_v = self.decoder_z1_z2(torch.cat((z2, ys), 1))

        # Sampling
        ql_m, ql_v, library = self.l_encoder(xs_)  # let's keep that ind. of y

        px_scale, px_rate, px_dropout = self.decoder(self.dispersion, z1, library, batch_index)

        reconst_loss = -log_zinb_positive(xs, px_rate, torch.exp(self.px_r), px_dropout)

        # KL Divergence
        kl_divergence_z2 = torch.sum(0.5 * (qz2_m ** 2 + qz2_v - torch.log(qz2_v + 1e-8) - 1), dim=1)
        kl_divergence_z1 = torch.sum(0.5 * (((qz1_m - pz1_m) ** 2 + qz1_v) / pz1_v - torch.log(qz1_v + 1e-8)
                                            + torch.log(pz1_v + 1e-8) - 1), dim=1)
        kl_divergence_l = torch.sum(0.5 * (((ql_m - local_l_mean) ** 2) / local_l_var + ql_v / local_l_var
                                           + torch.log(local_l_var + 1e-8) - torch.log(ql_v + 1e-8) - 1), dim=1)

        kl_divergence = (kl_divergence_z2 + kl_divergence_z1 + kl_divergence_l)

        if is_labelled:
            return reconst_loss, kl_divergence

        reconst_loss = reconst_loss.view(self.n_labels, -1)
        kl_divergence = kl_divergence.view(self.n_labels, -1)

        proba = self.classifier(z1_)

        reconst_loss = (reconst_loss.t() * proba).sum(dim=1)
        kl_divergence = (kl_divergence.t() * proba).sum(dim=1)
        y_prior = self.y_prior.type(proba.type())
        kl_divergence += torch.sum(torch.mul(proba, torch.log(y_prior) - torch.log(proba + 1e-8)), dim=-1)

        return reconst_loss, kl_divergence
