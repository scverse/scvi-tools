import time

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from scvi.external.harreman._utils import _resolve_device
from scvi.external.harreman.preprocessing.anndata import counts_from_anndata
from scvi.external.harreman.tools.knn import make_weights_non_redundant

from .local_autocorrelation import standardize_counts


def compute_local_correlation(
    adata: AnnData,
    genes: list | None = None,
    permutation_test: bool | None = False,
    M: int | None = 1000,
    seed: int | None = 42,
    check_analytic_null: bool | None = False,
    device: torch.device | str = "auto",
    verbose: bool | None = False,
):
    """
    Computes pairwise local correlation between selected genes using a spatial weight matrix.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData). Required keys in `adata.uns`:
        - 'gene_autocorrelation_results' (if `genes` is None)
        - 'layer_key': layer in `adata` to extract counts from
        - 'model': statistical model to use for centering (e.g., 'DANB', 'normal')
        - 'umi_counts': per-cell UMI counts
        - optionally 'sample_key': key in `adata.obs` to use for per-sample normalization
    genes : list, optional
        List of genes to include in the correlation analysis. If None, selects genes with
        FDR < 0.05 from `adata.uns['gene_autocorrelation_results']`, ordered by Z-score.
    permutation_test : bool, optional (default: False)
        Whether to compute an empirical p-value and null distribution by permuting the data.
    M : int, optional (default: 1000)
        Number of permutations to perform if `permutation_test` is True.
    seed : int, optional (default: 42)
        Random seed for permutation reproducibility.
    check_analytic_null : bool, optional (default: False)
        Whether to compute an analytic null distribution for the local correlation scores.
    device : torch.device, optional
        PyTorch device to run computations on. Defaults to CUDA if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        Results are stored in the following keys in `adata.uns`: `lcs`, `lc_zs`,
        `lc_z_pvals`, and `lc_z_FDR`.
    """
    start = time.time()

    if genes is None:
        gene_autocorrelation_results = adata.uns["gene_autocorrelation_results"]
        genes = (
            gene_autocorrelation_results.loc[gene_autocorrelation_results.Z_FDR < 0.05]
            .sort_values("Z", ascending=False)
            .index
        )

    if verbose:
        print(f"Computing pair-wise local correlation on {len(genes)} features...")

    layer_key = adata.uns["layer_key"]
    model = adata.uns["model"]
    sample_specific = "sample_key" in adata.uns.keys()

    # Load counts
    counts = counts_from_anndata(adata[:, genes], layer_key, dense=True)

    # UMI counts
    num_umi = adata.uns["umi_counts"]

    device = _resolve_device(device)

    # Convert to tensors
    num_umi = torch.tensor(adata.uns["umi_counts"], dtype=torch.float64, device=device)
    counts = torch.tensor(counts, dtype=torch.float64, device=device)

    # Center values
    counts = standardize_counts(adata, counts, model, num_umi, sample_specific)

    # Compute weights
    weights = make_weights_non_redundant(adata.obsp["weights"]).tocoo()
    weights = torch.sparse_coo_tensor(
        torch.tensor(np.vstack((weights.row, weights.col)), dtype=torch.long, device=device),
        torch.tensor(weights.data, dtype=torch.float64, device=device),
        torch.Size(weights.shape),
        device=device,
    )

    # Compute node degree
    row_degrees = torch.sparse.sum(weights, dim=1).to_dense()
    col_degrees = torch.sparse.sum(weights, dim=0).to_dense()
    D = row_degrees + col_degrees

    # Pairwise correlation
    WXt = torch.sparse.mm(weights, counts.T)  # (cells x genes)
    WtXt = torch.sparse.mm(weights.transpose(0, 1), counts.T)  # (cells x genes)
    lcs = torch.matmul(counts, WXt) + torch.matmul(counts, WtXt)

    # Compute second moments of H
    eg2s = (WXt + WtXt).pow(2).sum(dim=0)

    # Results
    results = compute_pairwise_correlation_results(
        counts=counts,
        weights=weights,
        lcs=lcs,
        eg2s=eg2s,
        genes=genes,
        D=D,
        M=M,
        permutation_test=permutation_test,
        seed=seed,
        check_analytic_null=check_analytic_null,
        device=device,
    )

    # Save results
    for key, value in results.items():
        adata.uns[key] = value

    if verbose:
        print(
            f"Pair-wise local correlation results are stored in adata.uns with the "
            f"following keys: {list(results.keys())}"
        )

        print(
            "Finished computing pair-wise local correlation in %.3f seconds"
            % (time.time() - start)
        )

    return


def compute_pairwise_correlation_results(
    counts, weights, lcs, eg2s, genes, D, M, permutation_test, seed, check_analytic_null, device
):
    """Compute pairwise local correlation results dict with Z-scores, p-values, and FDR."""
    results = {}

    lc_zs = compute_cor_Z_scores_torch(lcs, eg2s)

    lc_z_pvals = norm.sf(lc_zs.cpu().numpy())
    lc_z_FDR = multipletests(lc_z_pvals.flatten(), method="fdr_bh")[1].reshape(lc_z_pvals.shape)

    gene_index = pd.Index(genes)

    if permutation_test:
        n_genes, n_cells = counts.shape
        perm_array = torch.empty((n_genes, n_genes, M), dtype=torch.float16, device=device)
        if check_analytic_null:
            lc_zs_perm_array = torch.empty_like(perm_array)
            lc_pvals_perm_array = torch.empty_like(perm_array)

        torch.manual_seed(seed)
        for i in tqdm(range(M), desc="Permutation test"):
            idx = torch.randperm(n_cells, device=device)

            WXt_perm = torch.sparse.mm(weights, counts[:, idx].T)
            WtXt_perm = torch.sparse.mm(weights.transpose(0, 1), counts[:, idx].T)
            lcs_perm = torch.matmul(counts, WXt_perm) + torch.matmul(counts, WtXt_perm)

            perm_array[:, :, i] = lcs_perm.half()

            if check_analytic_null:
                lc_zs_perm = compute_cor_Z_scores_torch(lcs_perm, eg2s)
                lc_zs_perm_array[:, :, i] = lc_zs_perm.half()
                lc_pvals_perm_array[:, :, i] = torch.tensor(
                    norm.sf(lc_zs_perm.cpu().numpy()), device=device
                ).half()

        x = (perm_array > lcs.unsqueeze(-1)).sum(dim=2)
        lc_perm_pvals = (x + 1).float() / (M + 1)

        lc_perm_pvals_ab = torch.tril(lc_perm_pvals, diagonal=-1)
        lc_perm_pvals_ba = torch.tril(lc_perm_pvals.transpose(0, 1), diagonal=-1)
        lc_perm_pvals_sym = torch.where(
            lc_perm_pvals_ab > lc_perm_pvals_ba, lc_perm_pvals_ab, lc_perm_pvals_ba
        )
        i_upper = torch.triu_indices(n_genes, n_genes, offset=1)
        lc_perm_pvals_sym[i_upper[0], i_upper[1]] = lc_perm_pvals_sym[i_upper[1], i_upper[0]]

        results["lc_perm_pvals"] = pd.DataFrame(
            lc_perm_pvals.cpu().numpy(), index=gene_index, columns=gene_index
        )
        results["lc_perm_pvals_sym"] = pd.DataFrame(
            lc_perm_pvals_sym.cpu().numpy(), index=gene_index, columns=gene_index
        )

        if check_analytic_null:
            results["analytic_null_lc_zs_perm"] = lc_zs_perm_array.cpu().numpy()
            results["analytic_null_lc_pvals_perm"] = lc_pvals_perm_array.cpu().numpy()

    gene_maxs = 0.5 * torch.sum((counts**2) * D[None, :], dim=1)  # shape (n_genes,)
    lc_maxs = (gene_maxs[:, None] + gene_maxs[None, :]) / 2
    lcs = lcs / lc_maxs

    results["lcs"] = pd.DataFrame(lcs.cpu().numpy(), index=gene_index, columns=gene_index)
    results["lc_zs"] = pd.DataFrame(lc_zs.cpu().numpy(), index=gene_index, columns=gene_index)
    results["lc_z_pvals"] = pd.DataFrame(lc_z_pvals, index=gene_index, columns=gene_index)
    results["lc_z_FDR"] = pd.DataFrame(lc_z_FDR, index=gene_index, columns=gene_index)

    return results


def compute_cor_Z_scores_torch(lc, eg2s):
    """Compute symmetrized pairwise correlation Z-scores using PyTorch tensors."""
    EG = 0.0
    stdG = (eg2s - EG**2) ** 0.5

    Z = (lc - EG) / stdG

    Z_ab = torch.tril(Z, diagonal=-1)
    Z_ba = torch.tril(Z.T, diagonal=-1)

    Z = torch.where(Z_ab.abs() < Z_ba.abs(), Z_ab, Z_ba)
    Z = Z + Z.T

    return Z
