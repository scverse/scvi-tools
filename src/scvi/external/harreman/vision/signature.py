"""VISION-style gene signature scoring.

Adapted from visionpy (MIT licence, https://github.com/yoseflab/visionpy).
"""

from __future__ import annotations

import logging
import time
from re import compile, match
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix, issparse

from ._normalization import _get_mean_var, get_normalized_copy_sparse

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

    from ._normalization import NormData

logger = logging.getLogger(__name__)

_DOWN_SIG_KEY = "DN"
_UP_SIG_KEY = "UP"


# ---------------------------------------------------------------------------
# GPU-accelerated batch scoring
# ---------------------------------------------------------------------------


def _resolve_torch_device(device: str) -> str | None:
    """Return a torch device string if GPU is available, else ``None``."""
    try:
        import torch
    except ImportError:
        return None

    if device == "auto":
        if torch.cuda.is_available():
            return "cuda"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
        return None

    return device


def _sig_to_dense_f32(sig_matrix) -> np.ndarray:
    if isinstance(sig_matrix, pd.DataFrame):
        return sig_matrix.values.astype(np.float32)
    if sparse.issparse(sig_matrix):
        return sig_matrix.toarray().astype(np.float32)
    return np.asarray(sig_matrix, dtype=np.float32)


def _sig_to_scipy_csc(sig_matrix) -> sparse.csc_matrix:
    if isinstance(sig_matrix, pd.DataFrame):
        return sparse.csc_matrix(sig_matrix.values)
    if sparse.issparse(sig_matrix):
        return sig_matrix.tocsc()
    return sparse.csc_matrix(np.asarray(sig_matrix, dtype=float))


def _batch_sig_eval_norm_scipy(norm_data: NormData, sig_matrix, batch_size: int) -> np.ndarray:
    X = norm_data.data
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)

    gof = norm_data.gene_offsets
    gsf = norm_data.gene_scale_factors
    cof = norm_data.cell_offsets
    csf = norm_data.cell_scale_factors

    S_all = _sig_to_scipy_csc(sig_matrix)
    n_sigs = S_all.shape[1]
    chunks = []

    for start in range(0, n_sigs, batch_size):
        S = S_all[:, start : start + batch_size]
        S_gsf = S.multiply(gsf[:, None])
        t1 = (X @ S_gsf).toarray()
        t2 = np.asarray(S_gsf.multiply(gof[:, None]).sum(axis=0)).ravel()
        sig_sums = np.asarray(S.sum(axis=0)).ravel()
        t3 = np.outer(cof, sig_sums)
        result = (t1 + t2[None, :] + t3) * csf[:, None]
        denom = np.asarray(np.abs(S).sum(axis=0)).ravel()
        denom = np.where(denom > 0, denom, 1.0)
        chunks.append(result / denom[None, :])

    return np.concatenate(chunks, axis=1)


def _batch_sig_eval_norm_torch(
    norm_data: NormData, sig_matrix, device: str, batch_size: int
) -> np.ndarray:
    import torch

    dev = torch.device(device)
    use_cuda = dev.type == "cuda"

    X = norm_data.data
    if use_cuda:
        Xc = (X if sparse.issparse(X) else sparse.csr_matrix(X)).tocsr().astype(np.float32)
        X_pt = torch.sparse_csr_tensor(
            torch.from_numpy(Xc.indptr.copy().astype(np.int64)).to(dev),
            torch.from_numpy(Xc.indices.copy().astype(np.int64)).to(dev),
            torch.from_numpy(Xc.data.copy()).to(dev),
            size=tuple(Xc.shape),
        )
    else:
        X_dense = (X.toarray() if sparse.issparse(X) else np.asarray(X)).astype(np.float32)
        X_pt = torch.as_tensor(X_dense, device=dev)

    gof = torch.as_tensor(norm_data.gene_offsets.astype(np.float32), device=dev)
    gsf = torch.as_tensor(norm_data.gene_scale_factors.astype(np.float32), device=dev)
    cof = torch.as_tensor(norm_data.cell_offsets.astype(np.float32), device=dev)
    csf = torch.as_tensor(norm_data.cell_scale_factors.astype(np.float32), device=dev)

    S_full = torch.as_tensor(_sig_to_dense_f32(sig_matrix), device=dev)
    n_sigs = S_full.shape[1]

    chunks = []
    for start in range(0, n_sigs, batch_size):
        S = S_full[:, start : start + batch_size]
        S_gsf = S * gsf.unsqueeze(1)
        t1 = torch.sparse.mm(X_pt, S_gsf) if use_cuda else X_pt @ S_gsf
        t2 = (S_gsf * gof.unsqueeze(1)).sum(dim=0)
        t3 = torch.outer(cof, S.sum(dim=0))
        result = (t1 + t2.unsqueeze(0) + t3) * csf.unsqueeze(1)
        denom = S.abs().sum(dim=0).clamp(min=1e-10)
        chunks.append((result / denom.unsqueeze(0)).cpu())

    return torch.cat(chunks, dim=1).numpy()


def _batch_sig_eval_norm(
    norm_data: NormData,
    sig_matrix: np.ndarray | sparse.spmatrix | pd.DataFrame,
    device: str = "auto",
    batch_size: int = 1200,
) -> np.ndarray:
    """Score signatures via NormData without materialising the dense expression matrix."""
    torch_device = _resolve_torch_device(device)

    if torch_device is not None:
        try:
            return _batch_sig_eval_norm_torch(norm_data, sig_matrix, torch_device, batch_size)
        except Exception as exc:  # noqa: BLE001
            logger.warning("PyTorch path failed (%s); falling back to scipy sparse.", exc)

    return _batch_sig_eval_norm_scipy(norm_data, sig_matrix, batch_size)


# ---------------------------------------------------------------------------
# GMT / dict parsing
# ---------------------------------------------------------------------------


def _read_gmt(
    gmt_file: str,
    up_suffix: str = "(_up|_plus)",
    down_suffix: str = "(_down|_dn|_minus|_mius)",
    verbose: bool = False,
) -> dict[str, dict[str, list[str]]]:
    """Parse a GMT file and return a signed gene-set dict."""
    with open(gmt_file) as f:
        text_gmt = f.read().split("\n")

    pattern_down = compile(r"(\S+)" + down_suffix + "$")
    pattern_up = compile(r"(\S+)" + up_suffix + "$")
    signed_sign: dict[str, dict[str, list[str]]] = {}

    for line in text_gmt:
        temp_split = line.split("\t")
        sig_full_name = temp_split[0]
        if len(temp_split) < 3 or not sig_full_name:
            if verbose:
                logger.debug("Skipping: %s", sig_full_name)
            continue
        z = match(pattern_down, sig_full_name.lower())
        if z:
            sig_name, direction = z.groups()[0], _DOWN_SIG_KEY
        else:
            z = match(pattern_up, sig_full_name.lower())
            sig_name = z.groups()[0] if z else sig_full_name
            direction = _UP_SIG_KEY
        gene_list = [x for x in temp_split[2:] if x]
        if sig_name in signed_sign:
            signed_sign[sig_name][direction] = gene_list
        else:
            signed_sign[sig_name] = {direction: gene_list}

    return signed_sign


def _read_dict(
    d: dict,
    up_suffix: str = "(_up|_plus)",
    down_suffix: str = "(_down|_dn|_minus|_mius)",
) -> dict[str, dict[str, list[str]]]:
    """Parse a Python dict and return a signed gene-set dict."""
    pattern_down = compile(r"(\S+)" + down_suffix + "$")
    pattern_up = compile(r"(\S+)" + up_suffix + "$")
    signed_sign: dict[str, dict[str, list[str]]] = {}

    for sig_full_name, gene_list in d.items():
        z = match(pattern_down, sig_full_name.lower())
        if z:
            sig_name, direction = z.groups()[0], _DOWN_SIG_KEY
        else:
            z = match(pattern_up, sig_full_name.lower())
            sig_name = z.groups()[0] if z else sig_full_name
            direction = _UP_SIG_KEY
        if sig_name in signed_sign:
            signed_sign[sig_name][direction] = gene_list
        else:
            signed_sign[sig_name] = {direction: gene_list}

    return signed_sign


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_signatures(
    adata: AnnData,
    gmt_files: Sequence[str] | None = None,
    dicts: Sequence[dict] | None = None,
    use_raw: bool = False,
    min_signature_genes: int = 5,
    sig_gene_threshold: float = 0.001,
) -> None:
    r"""Load gene signatures into ``adata.varm["signatures"]``.

    Parses gene-set definitions from GMT files or Python dicts and builds a
    gene x signature weight matrix (+1 / -1 / 0) stored in
    ``adata.varm["signatures"]``.  Signatures with fewer than
    ``min_signature_genes`` genes matching the dataset are dropped.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.  Gene names are matched case-insensitively
        against the signature gene lists.
    gmt_files : list of str, optional
        Paths to one or more GMT files.  Each line must be tab-separated
        (name, description, gene1, gene2, ...).  Directional suffixes
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
        (n_genes x n_sigs) :class:`pandas.DataFrame` of +/-1 / 0 values.
    """
    if dicts is None and gmt_files is None:
        raise ValueError("Provide at least one of gmt_files or dicts.")

    files: list = []
    if dicts is not None:
        files.extend(dicts)
    if gmt_files is not None:
        files.extend(gmt_files)

    sig_dict: dict[str, dict[str, list[str]]] = {}
    for item in files:
        if isinstance(item, dict):
            sig_dict.update(_read_dict(item))
        else:
            sig_dict.update(_read_gmt(item))

    index = adata.raw.var.index if use_raw else adata.var_names
    columns = list(sig_dict.keys())
    sig_df = pd.DataFrame(0.0, index=index, columns=columns)
    sig_df.index = sig_df.index.str.lower()

    for sig_name, genes_up_down in sig_dict.items():
        for key, value in ((_UP_SIG_KEY, 1.0), (_DOWN_SIG_KEY, -1.0)):
            genes = genes_up_down.get(key)
            if genes is not None:
                genes_idx = pd.Index(genes).str.lower().intersection(sig_df.index)
                sig_df.loc[genes_idx, sig_name] = value

    sig_df.index = index

    if sig_gene_threshold > 0:
        X = adata.raw.X if use_raw else adata.X
        expressed_frac = (
            np.asarray((X > 0).mean(axis=0)).ravel()
            if issparse(X)
            else (np.asarray(X) > 0).mean(axis=0).ravel()
        )
        sig_df.iloc[expressed_frac < sig_gene_threshold, :] = 0.0

    sig_df = sig_df.loc[:, (sig_df != 0).any(axis=0)]
    if min_signature_genes > 0:
        n_matched = (sig_df != 0).sum(axis=0)
        sig_df = sig_df.loc[:, n_matched >= min_signature_genes]
        if sig_df.shape[1] == 0:
            raise ValueError(
                f"No signatures retained after requiring >= {min_signature_genes} "
                "matching genes. Lower min_signature_genes or use signatures with "
                "more gene overlap with this dataset."
            )

    adata.varm["signatures"] = sig_df


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
    n_cells x n_genes matrix: a :class:`~._normalization.NormData` container
    is built and scored in GPU-accelerated batches (CUDA -> MPS -> CPU
    fallback).  For dense inputs the standard VISION z-score formula is used.

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
        Key in ``adata.varm`` for the (n_genes x n_sigs) signature weight
        matrix (+1 upregulated, -1 downregulated, 0 absent).
    signature_names_uns_key : str or None, default None
        Key in ``adata.uns`` for custom signature column names.  When
        ``None``, names are read from the DataFrame columns (if applicable)
        or auto-generated as ``signature_0``, ``signature_1``, etc.
    sig_norm_method : str, default "znorm_columns"
        Normalisation method for the sparse scoring path.  Options:
        ``"znorm_columns"`` (z-normalise each cell across genes),
        ``"znorm_rows"`` (z-normalise each gene across cells),
        ``"znorm_rows_then_columns"``, ``"none"``.
    device : str, default "auto"
        Compute device: ``"auto"`` selects CUDA -> MPS -> CPU automatically.
    batch_size : int, default 1200
        Number of signatures scored per device batch.

    Returns
    -------
    None
        Results are stored in-place on ``adata``:

        - ``obsm["vision_signatures"]``: (n_cells x n_sigs) DataFrame of
          per-cell signature scores.
        - ``uns["norm_data_key"]``: records the expression key used.
        - ``uns["signature_varm_key"]``: records the varm key used.
    """
    t0 = time.time()
    logger.info("Computing signatures for each cell...")

    adata.uns["norm_data_key"] = norm_data_key
    adata.uns["signature_varm_key"] = signature_varm_key

    use_raw = norm_data_key == "use_raw"
    if norm_data_key is None:
        gene_expr = adata.X
    elif use_raw:
        gene_expr = adata.raw.X
    else:
        gene_expr = adata.layers[norm_data_key]

    sig_mat_raw = (
        adata.varm[signature_varm_key] if not use_raw else adata.raw.varm[signature_varm_key]
    )

    if signature_names_uns_key is not None:
        cols = adata.uns[signature_names_uns_key]
    elif isinstance(sig_mat_raw, pd.DataFrame):
        cols = sig_mat_raw.columns.tolist()
    else:
        cols = [f"signature_{i}" for i in range(sig_mat_raw.shape[1])]

    if issparse(gene_expr):
        norm_data = get_normalized_copy_sparse(gene_expr, method=sig_norm_method)
        scores = _batch_sig_eval_norm(norm_data, sig_mat_raw, device=device, batch_size=batch_size)
        sig_df = pd.DataFrame(data=scores, columns=cols, index=adata.obs_names)
    else:
        gene_expr = np.asarray(gene_expr, dtype=float)
        sig_matrix = csr_matrix(
            sig_mat_raw.to_numpy()
            if isinstance(sig_mat_raw, pd.DataFrame)
            else (sig_mat_raw.toarray() if issparse(sig_mat_raw) else sig_mat_raw)
        )
        cell_sig = gene_expr @ sig_matrix
        cell_sig = cell_sig.toarray() if issparse(cell_sig) else cell_sig
        sig_df = pd.DataFrame(data=cell_sig, columns=cols, index=adata.obs_names)
        mean, var = _get_mean_var(gene_expr, axis=1)
        n = np.asarray((sig_matrix > 0).sum(0))
        m = np.asarray((sig_matrix < 0).sum(0))
        denom = np.where((n + m) > 0, (n + m), 1.0)
        sig_df = sig_df / denom
        sig_mean = np.outer(mean, np.where(denom > 0, (n - m) / denom, 0.0))
        sig_std = np.sqrt(np.outer(var, np.where(denom > 0, 1.0 / denom, 1.0)))
        sig_std[sig_std == 0] = 1.0
        sig_df = (sig_df - sig_mean) / sig_std

    adata.obsm["vision_signatures"] = sig_df
    logger.info("Signatures computed in %.3f s.", time.time() - t0)
    return sig_df
