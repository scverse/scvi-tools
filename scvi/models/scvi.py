import logging
from anndata import AnnData

from scvi._compat import Literal
from scvi.core.models import VAE
from scvi.models._base import (
    BaseModelClass,
    RNASeqMixin,
    VAEMixin,
)

from scvi.core.trainers import UnsupervisedTrainer
from scvi.core.posteriors import Posterior

logger = logging.getLogger(__name__)


class SCVI(RNASeqMixin, VAEMixin, BaseModelClass):
    """
    single-cell Variational Inference [Lopez18]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for neural networks
    dispersion
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.dataset.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.models.SCVI(adata)
    >>> vae.train(n_epochs=400)
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()

    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        use_cuda: bool = True,
        **model_kwargs,
    ):
        super(SCVI, self).__init__(adata, use_cuda=use_cuda)
        self.model = VAE(
            n_input=self.summary_stats["n_genes"],
            n_batch=self.summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            reconstruction_loss=gene_likelihood,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self._model_summary_string = (
            "SCVI Model with following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            gene_likelihood,
            latent_distribution,
        )

    @property
    def _trainer_class(self):
        return UnsupervisedTrainer

    @property
    def _posterior_class(self):
        return Posterior
