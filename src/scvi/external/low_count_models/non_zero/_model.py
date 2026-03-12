"""Non-zero SCVI model."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scvi.model._scvi import SCVI

from ._module import NonZeroVAE

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


class NonZeroSCVI(SCVI):
    """SCVI variant that only computes reconstruction loss on non-zero observations.

    This model addresses the observation that zeros dominate training gradients in
    standard SCVI. With ~92-95% zeros in typical scRNA-seq data, the reconstruction
    loss is dominated by zero predictions, causing the model to hedge toward lower
    predictions.

    By masking zeros from the reconstruction loss, this model focuses on accurately
    predicting non-zero expression values, which may improve:

    - NLL/gene (nz): Lower (better non-zero prediction)
    - log2(mu/x) bias: Closer to 0 (less systematic under-prediction)
    - |log2(mu/x)|: Lower (better overall accuracy)

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
    normalize_by_nonzero
        If ``True``, normalize reconstruction loss by the number of non-zero
        observations per cell. This maintains a per-gene scale comparable to
        standard SCVI.
    **kwargs
        Additional keyword arguments passed to :class:`~scvi.model.SCVI`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.NonZeroSCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.external.NonZeroSCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()

    See Also
    --------
    :class:`~scvi.model.SCVI`
    :class:`~scvi.external.low_count_models.non_zero.NonZeroVAE`
    """

    _module_cls = NonZeroVAE

    def __init__(
        self,
        adata: AnnData | None = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson", "normal"] = "zinb",
        use_observed_lib_size: bool = True,
        latent_distribution: Literal["normal", "ln"] = "normal",
        normalize_by_nonzero: bool = True,
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
            normalize_by_nonzero=normalize_by_nonzero,
            **kwargs,
        )
        self._model_summary_string = (
            f"NonZeroSCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, latent_distribution: {latent_distribution}, "
            f"normalize_by_nonzero: {normalize_by_nonzero}."
        )
