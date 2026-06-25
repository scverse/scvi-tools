"""Sparse-preserving normalisation for VISION signature scoring.

Adapted from visionpy (MIT licence, https://github.com/yoseflab/visionpy).
Mirrors R VISION's NormalizationMethods.R.  Only the methods required for
signature scoring are included here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

# ---------------------------------------------------------------------------
# Tiny helpers (inlined from visionpy.utils and visionpy.projections)
# ---------------------------------------------------------------------------


def _get_mean_var(X, axis=0):
    if sparse.issparse(X):
        mean = np.array(X.mean(axis=axis)).ravel()
        var = np.array(X.power(2).mean(axis=axis)).ravel() - mean**2
    else:
        mean = np.mean(X, axis=axis)
        var = np.var(X, axis=axis, ddof=1)
    return mean, var


def _log2p1(X: np.ndarray | sparse.spmatrix) -> np.ndarray | sparse.spmatrix:
    """Sparse-preserving log2(x + 1) transform."""
    if sparse.issparse(X):
        X = X.copy()
        X.data = np.log2(X.data + 1)
    else:
        X = np.log2(np.asarray(X, dtype=float) + 1)
    return X


# ---------------------------------------------------------------------------
# NormData — lazy sparse-preserving normalisation container
# ---------------------------------------------------------------------------


@dataclass
class NormData:
    """Sparse expression matrix with pre-computed normalisation parameters.

    Stores the log-transformed matrix together with per-gene and per-cell
    offset / scale-factor vectors so that the fully normalised value for
    element ``(i, j)`` can be computed on demand::

        step1 = (data[i, j] + gene_offsets[j]) * gene_scale_factors[j]
        step2 = (step1 + cell_offsets[i]) * cell_scale_factors[i]

    Mirrors R VISION's ``NormData`` S4 class.
    """

    data: np.ndarray | sparse.spmatrix
    gene_offsets: np.ndarray
    gene_scale_factors: np.ndarray
    cell_offsets: np.ndarray
    cell_scale_factors: np.ndarray

    def toarray(self) -> np.ndarray:
        """Materialise the fully normalised dense (n_cells × n_genes) matrix."""
        X = (
            self.data.toarray()
            if sparse.issparse(self.data)
            else np.asarray(self.data, dtype=float)
        )
        X = (X + self.gene_offsets) * self.gene_scale_factors
        X = (X + self.cell_offsets[:, None]) * self.cell_scale_factors[:, None]
        return X


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _cell_norm_params_after_gene_norm(
    X: np.ndarray | sparse.spmatrix,
    gene_offsets: np.ndarray,
    gene_scale_factors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Analytically derive per-cell stats after gene normalisation.

    Mirrors R VISION's ``.colNormHelper``.  Avoids forming the dense
    n_cells × n_genes intermediate matrix.
    """
    n_cells, n_genes = X.shape
    gof = gene_offsets
    gsf = gene_scale_factors

    if sparse.issparse(X):
        xgsf_sum = np.asarray(X.multiply(gsf).sum(axis=1)).ravel()
    else:
        xgsf_sum = (X * gsf).sum(axis=1)

    const = float(np.sum(gof * gsf))
    cell_means = (xgsf_sum + const) / n_genes
    cell_offsets = -cell_means

    gsf2 = gsf**2
    gof_gsf2 = gof * gsf2

    if sparse.issparse(X):
        t1 = np.asarray(X.power(2).multiply(gsf2).sum(axis=1)).ravel()
        t2 = 2.0 * np.asarray(X.multiply(gof_gsf2).sum(axis=1)).ravel()
    else:
        t1 = ((X**2) * gsf2).sum(axis=1)
        t2 = 2.0 * (X * gof_gsf2).sum(axis=1)

    t3 = float(np.sum((gof**2) * gsf2))
    sum_sq = t1 + t2 + t3
    cell_var = (sum_sq - n_genes * cell_means**2) / max(n_genes - 1, 1)

    with np.errstate(invalid="ignore", divide="ignore"):
        cell_scale_factors = np.where(cell_var > 0, cell_var ** (-0.5), 1.0)

    return cell_offsets, cell_scale_factors


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

_SPARSE_METHODS = {"none", "znorm_columns", "znorm_rows", "znorm_rows_then_columns"}


def get_normalized_copy_sparse(
    X: np.ndarray | sparse.spmatrix,
    method: str = "znorm_columns",
) -> NormData:
    """Log2-transform then build a sparse-preserving :class:`NormData` object.

    Mirrors R VISION's ``getNormalizedCopySparse``.

    Parameters
    ----------
    X : array-like of shape (n_cells, n_genes)
        Expression matrix, cells × genes.
    method : str
        One of ``"none"``, ``"znorm_columns"`` (z-normalise each cell across
        genes), ``"znorm_rows"`` (z-normalise each gene across cells), or
        ``"znorm_rows_then_columns"``.

    Returns
    -------
    NormData
        Lazy normalisation container.
    """
    if method not in _SPARSE_METHODS:
        raise ValueError(
            f"Method '{method}' is not supported in the sparse path. "
            f"Choose from: {sorted(_SPARSE_METHODS)}."
        )

    X = _log2p1(X)
    n_cells, n_genes = X.shape

    gene_offsets = np.zeros(n_genes, dtype=float)
    gene_scale_factors = np.ones(n_genes, dtype=float)
    cell_offsets = np.zeros(n_cells, dtype=float)
    cell_scale_factors = np.ones(n_cells, dtype=float)

    if method in ("znorm_rows", "znorm_rows_then_columns"):
        gene_means, gene_var = _get_mean_var(X, axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            gsfact = np.where(gene_var > 0, gene_var ** (-0.5), 1.0)
        gene_offsets = -gene_means
        gene_scale_factors = gsfact

    if method == "znorm_columns":
        cell_means, cell_var = _get_mean_var(X, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            csfact = np.where(cell_var > 0, cell_var ** (-0.5), 1.0)
        cell_offsets = -cell_means
        cell_scale_factors = csfact

    if method == "znorm_rows_then_columns":
        cell_offsets, cell_scale_factors = _cell_norm_params_after_gene_norm(
            X, gene_offsets, gene_scale_factors
        )

    return NormData(
        data=X,
        gene_offsets=gene_offsets,
        gene_scale_factors=gene_scale_factors,
        cell_offsets=cell_offsets,
        cell_scale_factors=cell_scale_factors,
    )
