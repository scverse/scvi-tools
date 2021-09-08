import logging
from typing import List, Optional

from anndata import AnnData

from scvi._compat import Literal
from scvi.data._anndata import _setup_anndata
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.module import VAE

from .base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin

logger = logging.getLogger(__name__)


class SCVI(
    RNASeqMixin, VAEMixin, ArchesMixin, UnsupervisedTrainingMixin, BaseModelClass
):
    """
    single-cell Variational Inference [Lopez18]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
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
        Keyword args for :class:`~scvi.module.VAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.data.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/user_guide/notebooks/api_overview`
    2. :doc:`/user_guide/notebooks/harmonization`
    3. :doc:`/user_guide/notebooks/scarches_scvi_tools`
    4. :doc:`/user_guide/notebooks/scvi_in_R`
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
        **model_kwargs,
    ):
        super(SCVI, self).__init__(adata)

        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )
        self.module = VAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_labels=self.summary_stats["n_labels"],
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self._model_summary_string = (
            "SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
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
        self.init_params_ = self._get_init_params(locals())

    @staticmethod
    def setup_anndata(
        adata: AnnData,
        batch_key: Optional[str] = None,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
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
        labels_key
            key in `adata.obs` for label information. Categories will automatically be converted into integer
            categories and saved to `adata.obs['_scvi_labels']`. If `None`, assigns the same label to all the data.
        layer
            if not `None`, uses this as the key in `adata.layers` for raw count data. This is currently used only if a :class:`~scvi.model.SCANVI` model is created from this model.
        categorical_covariate_keys
            keys in `adata.obs` that correspond to categorical data.
        continuous_covariate_keys
            keys in `adata.obs` that correspond to continuous data.
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
            labels_key=labels_key,
            layer=layer,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys,
            copy=copy,
        )
