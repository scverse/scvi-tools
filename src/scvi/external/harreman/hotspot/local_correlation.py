from typing import Optional
import time
import numpy as np
import pandas as pd
import sparse
import torch
from anndata import AnnData
from numba import jit, njit
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm

from . import models
from ..preprocessing.anndata import counts_from_anndata
from .local_autocorrelation import compute_local_cov_max, standardize_counts
from ..tools.knn import make_weights_non_redundant


def compute_local_correlation(
    adata: AnnData,
    genes: Optional[list] = None,
    permutation_test: Optional[bool] = False,
    M: Optional[int] = 1000,
    seed: Optional[int] = 42,
    check_analytic_null: Optional[bool] = False,
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: Optional[bool] = False,
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
        List of genes to include in the correlation analysis. If None, selects genes with FDR < 0.05
        from `adata.uns['gene_autocorrelation_results']`, ordered by Z-score.
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
        Results are stored in the following keys in `adata.uns`: `lcs`, `lc_zs`, `lc_z_pvals`, and `lc_z_FDR`.
    """

    start = time.time()

    if genes is None:
        gene_autocorrelation_results = adata.uns['gene_autocorrelation_results']
        genes = gene_autocorrelation_results.loc[gene_autocorrelation_results.Z_FDR < 0.05].sort_values('Z', ascending=False).index

    if verbose:
        print(f"Computing pair-wise local correlation on {len(genes)} features...")
    
    layer_key = adata.uns['layer_key']
    model = adata.uns['model']
    sample_specific = 'sample_key' in adata.uns.keys()
    
    # Load counts
    counts = counts_from_anndata(adata[:, genes], layer_key, dense=True)
    
    # UMI counts
    num_umi = adata.uns['umi_counts']
    
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
        device=device)
    
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
        print(f"Pair-wise local correlation results are stored in adata.uns with the following keys: {list(results.keys())}")
    
        print("Finished computing pair-wise local correlation in %.3f seconds" %(time.time()-start))

    return


def compute_pairwise_correlation_results(
    counts, weights, lcs, eg2s, genes, D, M, permutation_test, seed, check_analytic_null, device
):
    
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
                lc_pvals_perm_array[:, :, i] = torch.tensor(norm.sf(lc_zs_perm.cpu().numpy()), device=device).half()
        
        x = (perm_array > lcs.unsqueeze(-1)).sum(dim=2)
        lc_perm_pvals = (x + 1).float() / (M + 1)
        
        lc_perm_pvals_ab = torch.tril(lc_perm_pvals, diagonal=-1)
        lc_perm_pvals_ba = torch.tril(lc_perm_pvals.transpose(0, 1), diagonal=-1)
        lc_perm_pvals_sym = torch.where(lc_perm_pvals_ab > lc_perm_pvals_ba, lc_perm_pvals_ab, lc_perm_pvals_ba)
        i_upper = torch.triu_indices(n_genes, n_genes, offset=1)
        lc_perm_pvals_sym[i_upper[0], i_upper[1]] = lc_perm_pvals_sym[i_upper[1], i_upper[0]]
        
        results["lc_perm_pvals"] = pd.DataFrame(lc_perm_pvals.cpu().numpy(), index=gene_index, columns=gene_index)
        results["lc_perm_pvals_sym"] = pd.DataFrame(lc_perm_pvals_sym.cpu().numpy(), index=gene_index, columns=gene_index)
        
        if check_analytic_null:
            results["analytic_null_lc_zs_perm"] = lc_zs_perm_array.cpu().numpy()
            results["analytic_null_lc_pvals_perm"] = lc_pvals_perm_array.cpu().numpy()

    gene_maxs = 0.5 * torch.sum((counts ** 2) * D[None, :], dim=1)  # shape (n_genes,)
    lc_maxs = (gene_maxs[:, None] + gene_maxs[None, :]) / 2
    lcs = lcs / lc_maxs
    
    results["lcs"] = pd.DataFrame(lcs.cpu().numpy(), index=gene_index, columns=gene_index)
    results["lc_zs"] = pd.DataFrame(lc_zs.cpu().numpy(), index=gene_index, columns=gene_index)
    results["lc_z_pvals"] = pd.DataFrame(lc_z_pvals, index=gene_index, columns=gene_index)
    results["lc_z_FDR"] = pd.DataFrame(lc_z_FDR, index=gene_index, columns=gene_index)
    
    return results


@jit(nopython=True)
def conditional_eg2(x, neighbors, weights):
    """
    Computes EG2 for the conditional correlation
    """
    N = neighbors.shape[0]
    K = neighbors.shape[1]
  
    t1x = np.zeros(N)

    for i in range(N):
        K = len(neighbors[i][~np.isnan(neighbors[i])])
        for k in range(K):
            j = neighbors[i, k]
            j = int(j)

            wij = weights[i, j]
            if wij == 0:
                continue

            t1x[i] += wij*x[j]
            t1x[j] += wij*x[i]

    out_eg2 = (t1x**2).sum()

    return out_eg2


@jit(nopython=True)
def local_cov_pair(x, y, neighbors, weights):
    """Test statistic for local pair-wise autocorrelation"""
    out = 0

    for i in range(len(x)):
        xi = x[i]
        yi = y[i]
        if xi == 0 and yi == 0:
            continue
        K = len(neighbors[i][~np.isnan(neighbors[i])])
        for k in range(K):

            j = neighbors[i, k]
            j = int(j)

            w_ij = weights[i, j]

            xj = x[j]
            yj = y[j]

            out += w_ij*(xi*yj + yi*xj)/2

    return out


@jit(nopython=True)
def local_cov_pair_fast(counts, weights):
    """Test statistic for local pair-wise autocorrelation"""
    counts_t = counts.transpose()
    weights_t = weights.transpose()

    lc_1 = sparse.einsum('ik,kl,lj->ij', counts, weights, counts_t)
    lc_2 = sparse.einsum('ik,kl,lj->ij', counts, weights_t, counts_t)
    lc = lc_1 + lc_2

    # lc = counts @ weights @ counts.T + counts @ weights.T @ counts.T

    return lc


def create_centered_counts_torch(counts, model, num_umi, device, eps=1e-10):
    """
    Vectorized PyTorch version of centered counts transformation.

    Args:
        counts: torch.Tensor of shape (G, C)
        model: one of {'bernoulli', 'danb', 'normal', 'none'}
        num_umi: torch.Tensor of shape (C,) = total UMIs per cell
        device

    Returns:
        centered_counts: torch.Tensor of shape (G, C)
    """
    G, C = counts.shape

    if model == 'bernoulli':
        vals = (counts > 0).double()
        model_fn = models.bernoulli_model_linear_torch
    elif model == 'danb':
        vals = counts.double()
        model_fn = models.danb_model_torch
    elif model == 'normal':
        vals = counts.double()
        model_fn = models.normal_model_torch
    elif model == 'none':
        # Center only (no variance normalization)
        mu = vals.mean(dim=1, keepdim=True)
        return counts - mu
    else:
        raise ValueError(f"Invalid model type: {model}")

    # Output tensors
    mu_all = torch.zeros_like(counts, dtype=torch.double)
    var_all = torch.ones_like(counts, dtype=torch.double)

    for g in range(G):
        mu, var, _ = model_fn(counts[g], num_umi, device)
        mu_all[g] = mu
        var_all[g] = var

    # Prevent division by zero
    var_all[var_all == 0] = 1.0

    centered = (counts - mu_all) / var_all.sqrt()
    centered[centered == 0] = 0.0  # Optional, usually unnecessary

    return centered


def create_centered_counts(counts, model, num_umi):
    """
    Creates a matrix of centered/standardized counts given
    the selected statistical model
    """
    out = np.zeros_like(counts, dtype='double')

    for i in range(out.shape[0]):

        vals_x = counts[i]

        out_x = create_centered_counts_row(
            vals_x, model, num_umi)

        out[i] = out_x

    return out


def create_centered_counts_row(vals_x, model, num_umi):

    if model == 'bernoulli':
        vals_x = (vals_x > 0).astype('double')
        mu_x, var_x, x2_x = models.bernoulli_model(
            vals_x, num_umi)
    elif model == 'danb':
        mu_x, var_x, x2_x = models.danb_model(
            vals_x, num_umi)
    elif model == 'normal':
        mu_x, var_x, x2_x = models.normal_model(
            vals_x, num_umi)
    elif model == 'none':
        mu_x, var_x, x2_x = models.none_model(
            vals_x, num_umi)
    else:
        raise Exception(f"Invalid Model: {model}")

    var_x[var_x == 0] = 1
    out_x = (vals_x-mu_x)/(var_x**0.5)
    out_x[out_x == 0] = 0

    return out_x


@jit(nopython=True)
def _compute_hs_pairs_inner_centered_cond_sym(
    rowpair, counts, neighbors, weights, eg2s
):
    """
    This version assumes that the counts have already been modeled
    and centered
    """
    row_i, row_j = rowpair

    vals_x = counts[row_i]
    vals_y = counts[row_j]

    lc = local_cov_pair(vals_x, vals_y, neighbors, weights)*2

    # Compute xy
    EG, EG2 = 0, eg2s[row_i]

    stdG = (EG2 - EG ** 2) ** 0.5

    Zxy = (lc - EG) / stdG

    # Compute yx
    EG, EG2 = 0, eg2s[row_j]

    stdG = (EG2 - EG ** 2) ** 0.5

    Zyx = (lc - EG) / stdG

    if abs(Zxy) < abs(Zyx):
        Z = Zxy
    else:
        Z = Zyx

    return (lc, Z)


def compute_cor_Z_scores(
    lc, eg2s
):

    EG, EG2 = 0, eg2s
    stdG = (EG2 - EG ** 2) ** 0.5

    Z = (lc - EG) / stdG[:, np.newaxis]

    Z_ab = np.tril(Z, k=-1)
    Z_ba = np.tril(Z.T, k=-1)

    Z = np.where(np.abs(Z_ab) < np.abs(Z_ba), Z_ab, Z_ba)

    i_upper = np.triu_indices(Z.shape[0], k=1)
    Z[i_upper] = Z.T[i_upper]

    return Z


def compute_cor_Z_scores_torch(
    lc, eg2s
):
    
    EG = 0.0
    stdG = (eg2s - EG**2)**0.5

    Z = (lc - EG) / stdG
    
    Z_ab = torch.tril(Z, diagonal=-1)
    Z_ba = torch.tril(Z.T, diagonal=-1)
    
    Z = torch.where(Z_ab.abs() < Z_ba.abs(), Z_ab, Z_ba)
    Z = Z + Z.T

    return Z


@njit
def expand_pairs(pairs, vals, N):

    out = np.zeros((N, N))

    for i in range(len(pairs)):

        x = pairs[i, 0]
        y = pairs[i, 1]
        v = vals[i]

        out[x, y] = v
        out[y, x] = v

    return out


def compute_max_correlation(node_degrees, counts):
    """
    For a Genes x Cells count matrix, compute the maximal pair-wise correlation
    between any two genes
    """

    N_GENES = counts.shape[0]

    gene_maxs = np.zeros(N_GENES)
    for i in range(N_GENES):
        gene_maxs[i] = compute_local_cov_max(node_degrees, counts[i])

    result = gene_maxs.reshape((-1, 1)) + gene_maxs.reshape((1, -1))
    result = result / 2
    return result

