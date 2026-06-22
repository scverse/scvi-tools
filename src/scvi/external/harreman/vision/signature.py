from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def load_signatures(
    adata: AnnData,
    gmt_files: list[str] | None = None,
    dicts: list[dict] | None = None,
    use_raw: bool = False,
    min_signature_genes: int = 5,
    sig_gene_threshold: float = 0.001,
) -> None:
    """Load gene signatures into ``adata.varm["signatures"]``.

    Parses gene-set definitions from GMT files or Python dicts and builds a
    gene × signature weight matrix (+1 / −1 / 0) stored in
    ``adata.varm["signatures"]``.  Signatures with fewer than
    ``min_signature_genes`` genes matching the dataset are dropped.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.  Gene names are matched case-insensitively
        against the signature gene lists.
    gmt_files : list of str, optional
        Paths to one or more GMT files.  Each line must be tab-separated
        (name, description, gene1, gene2, …).  Directional suffixes
        ``_up`` / ``_plus`` (upregulated) and ``_down`` / ``_dn`` /
        ``_minus`` (downregulated) are recognised; a signature name without
        a suffix is treated as upregulated.
    dicts : list of dict, optional
        Signature definitions as Python dicts.  Keys should follow the same
        directional-suffix convention, values are lists of gene symbols::

            {"pathway_up": ["GENE1", "GENE2"], "pathway_dn": ["GENE3"]}

    use_raw : bool, default False
        Use ``adata.raw`` gene names and expression matrix for gene-presence
        filtering.
    min_signature_genes : int, default 5
        Minimum number of genes from a signature that must be present in the
        dataset for the signature to be retained.
    sig_gene_threshold : float, default 0.001
        Genes expressed in fewer than this fraction of cells are excluded from
        signatures before the minimum-gene filter is applied.

    Returns
    -------
    None
        ``adata.varm["signatures"]`` is created or overwritten in place with a
        (n_genes × n_sigs) :class:`pandas.DataFrame` of ±1 / 0 values.
    """
    from visionpy.signature import signatures_from_file

    signatures_from_file(
        adata,
        use_raw=use_raw,
        gmt_files=gmt_files,
        dicts=dicts,
        min_signature_genes=min_signature_genes,
        sig_gene_threshold=sig_gene_threshold,
    )


def compute_vision_signatures(
    adata: AnnData,
    norm_data_key: str | None = None,
    signature_varm_key: str = "signatures",
    signature_names_uns_key: str | None = None,
    sig_norm_method: str = "znorm_columns",
    device: str = "auto",
    batch_size: int = 1200,
) -> None:
    """Compute per-cell VISION signature scores.

    For sparse expression matrices the computation avoids densifying the
    n_cells × n_genes matrix: a
    :class:`~visionpy.normalization.NormData` container is built and scored
    in GPU-accelerated batches (CUDA → MPS → CPU fallback).  For dense
    inputs the standard VISION z-score formula is used.

    Calling this function is a prerequisite for
    :func:`~scvi.external.harreman.hotspot.modules.integrate_vision_hotspot_results`.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.  Must contain ``adata.varm[signature_varm_key]``
        (produced by :func:`load_signatures`).
    norm_data_key : str or None, default None
        Expression matrix to score against.  ``None`` uses ``adata.X``,
        ``"use_raw"`` uses ``adata.raw.X``, any other string is looked up in
        ``adata.layers``.
    signature_varm_key : str, default "signatures"
        Key in ``adata.varm`` for the (n_genes × n_sigs) signature weight
        matrix (+1 upregulated, −1 downregulated, 0 absent).
    signature_names_uns_key : str or None, default None
        Key in ``adata.uns`` for custom signature column names.  When
        ``None``, names are read from the DataFrame columns (if applicable)
        or auto-generated as ``signature_0``, ``signature_1``, etc.
    sig_norm_method : str, default "znorm_columns"
        Normalisation method for the sparse scoring path.  Options:
        ``"znorm_columns"`` (z-normalise each cell across genes),
        ``"znorm_rows"`` (z-normalise each gene across cells),
        ``"znorm_rows_then_columns"``, ``"rank_norm_columns"``, ``"none"``.
    device : str, default "auto"
        Compute device: ``"auto"`` selects CUDA → MPS → CPU automatically.
    batch_size : int, default 1200
        Number of signatures scored per device batch.

    Returns
    -------
    None
        Results are stored in-place on ``adata``:

        - ``obsm["vision_signatures"]``: (n_cells × n_sigs) DataFrame of
          per-cell signature scores.
        - ``uns["norm_data_key"]``: records the expression key used.
        - ``uns["signature_varm_key"]``: records the varm key used.
    """
    from visionpy.signature import compute_signatures_anndata

    compute_signatures_anndata(
        adata,
        norm_data_key=norm_data_key,
        signature_varm_key=signature_varm_key,
        signature_names_uns_key=signature_names_uns_key,
        sig_norm_method=sig_norm_method,
        device=device,
        batch_size=batch_size,
    )
