import logging
from typing import Optional

import pandas as pd
from anndata import AnnData

from scvi._compat import Literal
from scvi.data._anndata import _setup_anndata
from scvi.model._utils import _get_var_names_from_setup_anndata
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.module import LDVAE

from .base import BaseModelClass, RNASeqMixin, VAEMixin

logger = logging.getLogger(__name__)


class LinearSCVI(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """
    Linearly-decoded VAE [Svensson20]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.LinearSCVI.setup_anndata`.
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
    **model_kwargs
        Keyword args for :class:`~scvi.module.LDVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.LinearSCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.LinearSCVI(adata)
    >>> vae.train()
    >>> adata.var["loadings"] = vae.get_loadings()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/user_guide/notebooks/linear_decoder`
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
        **model_kwargs,
    ):
        super(LinearSCVI, self).__init__(adata)
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

    @staticmethod
    def setup_anndata(
        adata: AnnData,
        batch_key: Optional[str] = None,
        layer: Optional[str] = None,
        copy: bool = False,
    ) -> Optional[AnnData]:
        """
        Sets up the :class:`~anndata.AnnData` object for this model.

        A mapping will be created between data fields used by this model to their respective locations in adata.
        This method will also compute the log mean and log variance per batch for the library size prior.

        None of the data in adata are modified. Only adds fields to adata.

        Parameters
        ----------
        adata
            AnnData object containing raw counts. Rows represent cells, columns represent features.
        batch_key
            key in `adata.obs` for batch information. Categories will automatically be converted into integer
            categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch to all the data.
        layer
            if not `None`, uses this as the key in `adata.layers` for raw count data.
        copy
            if `True`, a copy of adata is returned.

        Returns
        -------
        If ``copy``,  will return :class:`~anndata.AnnData`.
        Adds the following fields to adata:

        .uns['_scvi']
            `scvi` setup dictionary
        .obs['_local_l_mean']
            per batch library size mean
        .obs['_local_l_var']
            per batch library size variance
        .obs['_scvi_labels']
            labels encoded as integers
        .obs['_scvi_batch']
            batch encoded as integers
        """
        return _setup_anndata(
            adata,
            batch_key=batch_key,
            layer=layer,
            copy=copy,
        )
