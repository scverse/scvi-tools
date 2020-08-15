import logging
import torch
import pandas as pd
from anndata import AnnData

from scvi._compat import Literal
from scvi.models._modules.vae import LDVAE
from scvi.models import SCVI

from scvi.inference.posterior import Posterior

logger = logging.getLogger(__name__)


class LinearSCVI(SCVI):
    """Linearly-decoded VAE [Svensson20]_

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder NN
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
    >>> vae = scvi.models.LinearSCVI(adata)
    >>> vae.train(n_epochs=400)
    >>> adata.var["loadings"] = vae.get_loadings()
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        use_cuda: bool = True,
        **model_kwargs,
    ):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"

        self.adata = adata
        summary_stats = adata.uns["scvi_summary_stats"]
        self.model = LDVAE(
            n_input=summary_stats["n_genes"],
            n_batch=summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            reconstruction_loss=gene_likelihood,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self.is_trained = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self.batch_size = 128
        self._posterior_class = Posterior

    def get_loadings(self) -> pd.DataFrame:
        """Extract per-gene weights in the linear decoder

        Shape is genes by dim(Z)
        """

        cols = ["Z_{}".format(i) for i in self.n_latent]
        loadings = pd.DataFrame(
            self.model.get_loadings(), index=self.adata.var_names, columns=cols
        )

        return loadings
