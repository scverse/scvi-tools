import torch

from scvi.models.vae import VAE
from scvi.modules import Classifier
from scvi.utils import enumerate_discrete, one_hot


# VAE model - for classification: VAEC
class VAEC(VAE):
    def __init__(self, n_input, n_labels=7, y_prior=None, n_layers=3, **kargs):
        super(VAEC, self).__init__(n_input=n_input + n_labels, n_layers=n_layers, **kargs)
        self.n_labels = n_labels

        self.y_prior = y_prior if y_prior is not None else (1 / self.n_labels) * torch.ones(self.n_labels)
        self.classifier = Classifier(n_input, self.n_hidden, self.n_labels, n_layers=n_layers,
                                     dropout_rate=self.dropout_rate)

    def classify(self, x):
        return self.classifier(x)

    def sample_from_posterior_z(self, x, y):
        x = torch.cat((x, one_hot(y, self.n_labels)), 1)
        return super(VAEC, self).sample_from_posterior_z(x)

    def sample_from_posterior_l(self, x, y):
        x = torch.cat((x, one_hot(y, self.n_labels)), 1)
        return super(VAEC, self).sample_from_posterior_l(x)

    def get_sample_scale(self, x, y, batch_index=None):
        z = self.sample_from_posterior_z(x, y)
        px = self.decoder.px_decoder_batch(z, batch_index)
        px_scale = self.decoder.px_scale_decoder(px)
        return px_scale

    def get_sample_rate(self, x, y, batch_index=None):
        z = self.sample_from_posterior_z(x, y)
        library = self.sample_from_posterior_l(x, y)
        px = self.decoder.px_decoder_batch(z, batch_index)
        return self.decoder.px_scale_decoder(px) * torch.exp(library)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        # Prepare for sampling
        xs, ys = (x, y)

        # Enumerate choices of label
        if not is_labelled:
            ys = enumerate_discrete(xs, self.n_labels)
            xs = xs.repeat(self.n_labels, 1)
            if batch_index is not None:
                batch_index = batch_index.repeat(self.n_labels, 1)
            local_l_var = local_l_var.repeat(self.n_labels, 1)
            local_l_mean = local_l_mean.repeat(self.n_labels, 1)
        else:
            ys = one_hot(ys, self.n_labels)

        reconst_loss, kl_divergence = super(VAEC, self).forward(torch.cat((xs, ys), 1), local_l_mean, local_l_var,
                                                                batch_index=batch_index)

        if is_labelled:
            return reconst_loss, kl_divergence

        reconst_loss = reconst_loss.view(self.n_labels, -1)
        kl_divergence = kl_divergence.view(self.n_labels, -1)
        if self.log_variational:
            x_ = torch.log(1 + x)  # Necessarily goes from z ! No not necessarily ! goes from x
        proba = self.classify(x_)  # Here it passes through z since that's inherited ! That's the magic
        reconst_loss = (reconst_loss.t() * proba).sum(dim=1)
        kl_divergence = (kl_divergence.t() * proba).sum(dim=1)
        y_prior = self.y_prior.type(proba.type())
        kl_divergence += torch.sum(torch.mul(proba, torch.log(y_prior) - torch.log(proba + 1e-8)), dim=-1)

        return reconst_loss, kl_divergence
