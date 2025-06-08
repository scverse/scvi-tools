from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.model._utils import _init_library_size
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.module import LDVAE
from scvi.utils import setup_anndata_dsp

from .base import BaseMinifiedModeModelClass, RNASeqMixin, VAEMixin

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class LinearSCVI(
    RNASeqMixin,
    VAEMixin,
    UnsupervisedTrainingMixin,
    BaseMinifiedModeModelClass,
):
    """Linearly-decoded VAE :cite:p:`Svensson20`.

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
    use_observed_lib_size
        If ``True``, use the observed library size for RNA as the scaling factor in the mean of the
        conditional distribution.
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

    1. :doc:`/tutorials/notebooks/scrna/linear_decoder`
    """

    _module_cls = LDVAE
    _LATENT_QZM_KEY = "ldvae_latent_qzm"
    _LATENT_QZV_KEY = "ldvae_latent_qzv"

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        use_observed_lib_size: bool = False,
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch
        library_log_means, library_log_vars = None, None
        if (
            self.minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR
            and not use_observed_lib_size
        ):
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            use_observed_lib_size=use_observed_lib_size,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"LinearSCVI Model with the following params: \nn_hidden: {n_hidden}, "
            f"n_latent: {n_latent}, n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"dispersion: {dispersion}, gene_likelihood: {gene_likelihood}, "
            f"latent_distribution: {latent_distribution}"
        )
        self.n_latent = n_latent
        self.init_params_ = self._get_init_params(locals())

    def get_loadings(self) -> pd.DataFrame:
        """Extract per-gene weights in the linear decoder.

        Shape is genes by `n_latent`.

        """
        cols = [f"Z_{i}" for i in range(self.n_latent)]
        var_names = self.adata.var_names
        loadings = pd.DataFrame(self.module.get_loadings(), index=var_names, columns=cols)

        return loadings

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        layer: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
