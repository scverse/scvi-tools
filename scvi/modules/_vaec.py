import torch
from torch.distributions import Categorical, Normal
from torch.distributions import kl_divergence as kl

from scvi.compose import DecoderSCVI, Encoder
from scvi.modules import Classifier
from scvi.modules._utils import broadcast_labels

from ._vae import VAE


class VAEC(VAE):
    r"""
    A semi-supervised Variational auto-encoder model - inspired from M2 model.

    Described in (https://arxiv.org/pdf/1406.5298.pdf)

    Parameters
    ----------
    n_input :
        Number of input genes
    n_batch :
        Number of batches
    n_labels :
        Number of labels
    n_hidden :
        Number of nodes per hidden layer
    n_latent :
        Dimensionality of the latent space
    n_layers :
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate :
        Dropout rate for neural networks
    dispersion :
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    log_variational :
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood :
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    y_prior :
        If None, initialized to uniform probability over cell types

    Examples
    --------
    >>> gene_dataset = CortexDataset()
    >>> vaec = VAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
    ... n_labels=gene_dataset.n_labels)

    >>> gene_dataset = SyntheticDataset(n_labels=3)
    >>> vaec = VAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
    ... n_labels=3, y_prior=torch.tensor([[0.1,0.5,0.4]]))

    """

    def __init__(
        self,
        n_input,
        n_batch,
        n_labels,
        n_hidden=128,
        n_latent=10,
        n_layers=1,
        dropout_rate=0.1,
        y_prior=None,
        dispersion="gene",
        log_variational=True,
        gene_likelihood="zinb",
    ):
        super().__init__(
            n_input,
            n_batch,
            n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            use_observed_lib_size=False,
        )

        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_cat_list=[n_batch, n_labels],
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
        )
        self.decoder = DecoderSCVI(
            n_latent,
            n_input,
            n_cat_list=[n_batch, n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
        )

        self.y_prior = torch.nn.Parameter(
            y_prior
            if y_prior is not None
            else (1 / n_labels) * torch.ones(1, n_labels),
            requires_grad=False,
        )

        self.classifier = Classifier(
            n_input, n_hidden, n_labels, n_layers=n_layers, dropout_rate=dropout_rate
        )

    def classify(self, x, batch_index=None):
        x = torch.log(1 + x)
        return self.classifier(x)

    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        is_labelled = False if y is None else True

        # Prepare for sampling
        x_ = torch.log(1 + x)
        ql_m, ql_v, library = self.l_encoder(x_, batch_index)

        # Enumerate choices of label
        ys, xs, library_s, batch_index_s = broadcast_labels(
            y, x, library, batch_index, n_broadcast=self.n_labels
        )

        # Sampling
        outputs = self.inference(xs, batch_index_s, ys)
        px_r = outputs["px_r"]
        px_rate = outputs["px_rate"]
        px_dropout = outputs["px_dropout"]
        qz_m = outputs["qz_m"]
        qz_v = outputs["qz_v"]
        reconst_loss = self.get_reconstruction_loss(xs, px_rate, px_r, px_dropout)

        # KL Divergence
        mean = torch.zeros_like(qz_m)
        scale = torch.ones_like(qz_v)

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(mean, scale)).sum(
            dim=1
        )
        kl_divergence_l = kl(
            Normal(ql_m, torch.sqrt(ql_v)),
            Normal(local_l_mean, torch.sqrt(local_l_var)),
        ).sum(dim=1)

        if is_labelled:
            return reconst_loss, kl_divergence_z + kl_divergence_l, 0.0

        reconst_loss = reconst_loss.view(self.n_labels, -1)

        probs = self.classifier(x_)
        reconst_loss = (reconst_loss.t() * probs).sum(dim=1)

        kl_divergence = (kl_divergence_z.view(self.n_labels, -1).t() * probs).sum(dim=1)
        kl_divergence += kl(
            Categorical(probs=probs),
            Categorical(probs=self.y_prior.repeat(probs.size(0), 1)),
        )
        kl_divergence += kl_divergence_l

        return reconst_loss, kl_divergence, 0.0
