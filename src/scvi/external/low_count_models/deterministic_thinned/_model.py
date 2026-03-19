"""Deterministic Thinned SCVI model for ablation study."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scvi.model._scvi import SCVI

from ._module import DeterministicThinnedVAE

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


class DeterministicThinnedSCVI(SCVI):
    """ThinnedSCVI ablation using encoder mean instead of sampling.

    This model extends SCVI by training on thinned versions of the input data
    (like ThinnedSCVI) but uses the encoder mean directly instead of sampling
    from the posterior distribution.

    This is an ablation to test whether:
    1. The stochastic noise injection from the reparameterization trick
       is important for robustness during VAE training.
    2. Exposure to thinned data alone (via data augmentation) is sufficient
       for robustness to low-UMI cells.

    The model still:
    - Trains on randomly thinned data (data augmentation from ThinnedVAE)
    - Computes KL divergence loss (still regularizes the encoder)
    - Uses the same architecture as standard SCVI

    The key difference is that z = encoder_mean instead of z ~ N(encoder_mean, encoder_var).

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
    min_library_size
        Minimum target library size for thinning. Default is 10.
        Thinned library sizes are sampled log-uniformly between this
        value and the observed library size.
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
        * ``'normal'`` - ``EXPERIMENTAL`` Normal distribution
    use_observed_lib_size
        If ``True``, use the observed library size for RNA as the scaling factor in the mean of the
        conditional distribution.
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    **kwargs
        Additional keyword arguments for
        :class:`~scvi.external.low_count_models.deterministic_thinned.DeterministicThinnedVAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.DeterministicThinnedSCVI.setup_anndata(adata, batch_key="batch")
    >>> model = scvi.external.DeterministicThinnedSCVI(adata, min_library_size=10.0)
    >>> model.train()
    >>> latent = model.get_latent_representation()

    See Also
    --------
    :class:`~scvi.model.SCVI`
    :class:`~scvi.external.ThinnedSCVI`
    :class:`~scvi.external.low_count_models.deterministic_thinned.DeterministicThinnedVAE`
    """

    _module_cls = DeterministicThinnedVAE

    def __init__(
        self,
        adata: AnnData | None = None,
        min_library_size: float = 10.0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson", "normal"] = "zinb",
        use_observed_lib_size: bool = True,
        latent_distribution: Literal["normal", "ln"] = "normal",
        **kwargs,
    ):
        super().__init__(
            adata=adata,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            use_observed_lib_size=use_observed_lib_size,
            latent_distribution=latent_distribution,
            min_library_size=min_library_size,
            **kwargs,
        )
        self._model_summary_string = (
            "DeterministicThinnedSCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, latent_distribution: {latent_distribution}, "
            f"min_library_size: {min_library_size}."
        )
