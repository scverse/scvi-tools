"""Joint Embedding SCVI model."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scvi.model._scvi import SCVI

from ._module import JointEmbeddingVAE

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


class JointEmbeddingSCVI(SCVI):
    """SCVI with joint embedding loss using binomial thinning and CCO.

    This model extends the standard SCVI with a cross-correlation objective (CCO)
    loss that encourages the embedding of a thinned view to match the embedding
    of the original data. This promotes robustness to count dropout/noise,
    particularly for low-UMI cells that would otherwise converge to a learned
    bias point in the encoder.

    Thinning probabilities are dynamically sampled per cell to produce target
    library sizes that are log-uniform between min_library_size and the observed
    library size. This matches realistic library size variation in single-cell data.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
    joint_embedding_weight
        Weight for the CCO loss. Default is 1.0.
    lambda_off_diag
        Off-diagonal penalty in CCO loss. Default is 0.01.
    min_library_size
        Minimum target library size for thinning. Default is 10.
        Thinned library sizes are sampled log-uniformly between this
        value and the observed library size.
    reconstruction_weight
        Weight for reconstruction loss. Default is 1.0.
        Set to 0.0 for pure self-supervised training with only CCO loss.
    variance_weight
        Weight for variance regularization loss (VICReg-style). Default is 0.0.
        Set to positive value (e.g., 1.0) to prevent dimension collapse
        in self-supervised training.
    use_joint_embedding
        Whether to use joint embedding loss. Default is True.
        Set to False to train as standard SCVI.
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
        :class:`~scvi.external.low_count_models.joint_embedding.JointEmbeddingVAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.JointEmbeddingSCVI.setup_anndata(adata, batch_key="batch")
    >>> model = scvi.external.JointEmbeddingSCVI(
    ...     adata,
    ...     joint_embedding_weight=1.0,
    ...     lambda_off_diag=0.01,
    ... )
    >>> model.train()
    >>> latent = model.get_latent_representation()

    See Also
    --------
    :class:`~scvi.model.SCVI`
    :class:`~scvi.external.ThinnedSCVI`
    :class:`~scvi.external.low_count_models.joint_embedding.JointEmbeddingVAE`
    """

    _module_cls = JointEmbeddingVAE

    def __init__(
        self,
        adata: AnnData | None = None,
        joint_embedding_weight: float = 1.0,
        lambda_off_diag: float = 0.01,
        min_library_size: float = 10.0,
        reconstruction_weight: float = 1.0,
        variance_weight: float = 0.0,
        use_joint_embedding: bool = True,
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
            joint_embedding_weight=joint_embedding_weight,
            lambda_off_diag=lambda_off_diag,
            min_library_size=min_library_size,
            reconstruction_weight=reconstruction_weight,
            variance_weight=variance_weight,
            use_joint_embedding=use_joint_embedding,
            **kwargs,
        )
        self._model_summary_string = (
            "JointEmbeddingSCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, latent_distribution: {latent_distribution}, "
            f"joint_embedding_weight: {joint_embedding_weight}, "
            f"lambda_off_diag: {lambda_off_diag}, min_library_size: {min_library_size}, "
            f"reconstruction_weight: {reconstruction_weight}, "
            f"variance_weight: {variance_weight}, "
            f"use_joint_embedding: {use_joint_embedding}."
        )
