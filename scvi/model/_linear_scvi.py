import logging

import pandas as pd
from anndata import AnnData

from scvi._compat import Literal
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import TrainingPlan
from scvi.model._utils import _get_var_names_from_setup_anndata
from scvi.modules import LDVAE

from .base import BaseModelClass, RNASeqMixin, VAEMixin

logger = logging.getLogger(__name__)


class LinearSCVI(RNASeqMixin, VAEMixin, BaseModelClass):
    """
    Linearly-decoded VAE [Svensson20]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder NN.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.modules.LDVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.data.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.LinearSCVI(adata)
    >>> vae.train()
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
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(LinearSCVI, self).__init__(adata, use_gpu=use_gpu)
        self.module = LDVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self._model_summary_string = (
            "LinearSCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
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
        self.n_latent = n_latent
        self.init_params_ = self._get_init_params(locals())

    @property
    def _plan_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    def get_loadings(self) -> pd.DataFrame:
        """
        Extract per-gene weights in the linear decoder.

        Shape is genes by `n_latent`.

        """
        cols = ["Z_{}".format(i) for i in range(self.n_latent)]
        var_names = _get_var_names_from_setup_anndata(self.adata)
        loadings = pd.DataFrame(
            self.module.get_loadings(), index=var_names, columns=cols
        )

        return loadings
