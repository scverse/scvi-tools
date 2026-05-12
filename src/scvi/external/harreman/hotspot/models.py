from numba import jit
import pandas as pd
import numpy as np
from numba import njit
import torch
from typing import Callable, Union


def danb_model(gene_counts, umi_counts):

    tj = gene_counts.sum()
    tis = umi_counts
    total = tis.sum()

    N = gene_counts.size

    min_size = 10**(-10)

    mu = tj*tis/total
    vv = (gene_counts - mu).var()*(N/(N-1))
    # vv = ((gene_counts - mu)**2).sum()
    my_rowvar = vv

    size = ((tj**2) / total) * ((tis**2).sum() / total) / ((N-1)*my_rowvar-tj)
    # size = ((tj**2) * ((tis/total)**2).sum()) / ((N-1)*my_rowvar-tj)

    if size < 0:    # Can't have negative dispersion
        size = 1e9

    if size < min_size and size >= 0:
        size = min_size

    var = mu*(1+mu/size)
    x2 = var+mu**2

    return mu, var, x2


def danb_model_torch(counts: torch.Tensor, umi_counts: torch.Tensor, eps: float = 1e-10):
    """
    Vectorized DANB model computation in PyTorch for a batch of genes.

    Args:
        counts: Tensor of shape [genes, cells], gene expression counts.
        umi_counts: Tensor of shape [cells], total UMI per cell.
        eps: Small constant to avoid division by zero or log(0).

    Returns:
        mu: Mean per gene per cell [genes, cells]
        var: Variance per gene per cell [genes, cells]
        x2: Second moment per gene per cell [genes, cells]
    """
    tj = counts.sum(dim=1, keepdim=True)              # [genes, 1]
    total = umi_counts.sum()                          # scalar
    N = counts.shape[1]                               # number of cells

    mu = tj * umi_counts / total                      # [genes, cells]
    diff = counts - mu                                # [genes, cells]

    # Unbiased sample variance (N / (N - 1))
    var_gene = (diff ** 2).mean(dim=1) * N / (N - 1)  # [genes]

    numerator = ((tj ** 2) / total).squeeze() * (umi_counts ** 2).sum() / total  # [genes]
    denominator = (N - 1) * var_gene - tj.squeeze()                               # [genes]
    size = numerator / (denominator + eps)                                       # [genes]

    # Clamp size for numerical stability
    size = torch.where(size < 0, torch.tensor(1e9, device=size.device), size)
    size = torch.clamp(size, min=eps)

    size = size.unsqueeze(1)                      # [genes, 1] for broadcasting
    var = mu * (1 + mu / size)                    # [genes, cells]
    x2 = var + mu**2                              # [genes, cells]

    return mu, var, x2


def ct_danb_model(gene_counts, umi_counts, cell_types):

    mu_ct = np.zeros(len(cell_types))
    var_ct = np.zeros(len(cell_types))
    x2_ct = np.zeros(len(cell_types))
    
    min_size = 10**(-10)
    
    for cell_type in np.unique(cell_types):

        gene_counts_ct = gene_counts[cell_types == cell_type]
        umi_counts_ct = umi_counts[cell_types == cell_type]
    
        tj = gene_counts_ct.sum()
        tis = umi_counts_ct
        total = tis.sum()

        N = gene_counts_ct.size
        
        mu = tj*tis/total
        vv = (gene_counts_ct - mu).var()*(N/(N-1)) if N>1 else (gene_counts_ct - mu).var()
        my_rowvar = vv

        size = ((tj**2) / total) * ((tis**2).sum() / total) / ((N-1)*my_rowvar-tj)

        if size < 0:    # Can't have negative dispersion
            size = 1e9

        if size < min_size and size >= 0:
            size = min_size

        var = mu*(1+mu/size)
        x2 = var+mu**2

        mu_ct[cell_types == cell_type] = mu
        var_ct[cell_types == cell_type] = var
        x2_ct[cell_types == cell_type] = x2

    return mu_ct, var_ct, x2_ct


N_BIN_TARGET = 30


@jit(nopython=True)
def find_gene_p(num_umi, D):
    """
    Finds gene_p such that sum of expected detects
    matches our data

    Performs a binary search on p in the space of log(p)
    """

    low = 1e-12
    high = 1

    if D == 0:
        return 0

    for ITER in range(40):

        attempt = (high*low)**0.5
        tot = 0

        for i in range(len(num_umi)):
            tot = tot + 1-(1-attempt)**num_umi[i]

        if abs(tot-D)/D < 1e-3:
            break

        if tot > D:
            high = attempt
        else:
            low = attempt

    return (high*low)**0.5


def bernoulli_model_scaled(gene_detects, umi_counts):

    D = gene_detects.sum()

    gene_p = find_gene_p(umi_counts, D)

    detect_p = 1-(1-gene_p)**umi_counts

    mu = detect_p
    var = detect_p * (1 - detect_p)
    x2 = detect_p

    return mu, var, x2


def true_params_scaled(gene_p, umi_counts):

    detect_p = 1-(1-gene_p/10000)**umi_counts

    mu = detect_p
    var = detect_p * (1 - detect_p)
    x2 = detect_p

    return mu, var, x2


def bernoulli_model_linear(gene_detects, umi_counts):

    # We modify the 0 UMI counts to 1e-10 to remove the NaN values from the qcut output.
    umi_counts[umi_counts == 0] = 1e-10

    umi_count_bins, bins = pd.qcut(
        np.log10(umi_counts), N_BIN_TARGET, labels=False, retbins=True,
        duplicates='drop'
    )
    bin_centers = np.array(
        [bins[i] / 2 + bins[i + 1] / 2 for i in range(len(bins) - 1)]
    )

    N_BIN = len(bin_centers)

    bin_detects = bin_gene_detection(gene_detects, umi_count_bins, N_BIN)

    lbin_detects = logit(bin_detects)

    X = np.ones((N_BIN, 2))
    X[:, 1] = bin_centers
    Y = lbin_detects

    b = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    detect_p = ilogit(b[0] + b[1] * np.log10(umi_counts))

    mu = detect_p
    var = detect_p * (1 - detect_p)
    x2 = detect_p

    return mu, var, x2


def bernoulli_model_linear_torch(gene_detects, umi_counts, n_bins=30, eps=1e-10):
    """
    gene_detects: [genes, cells] binary tensor
    umi_counts: [cells] tensor
    Returns: mu, var, x2 = each of shape [genes, cells]
    """
    device = gene_detects.device
    log_umi = torch.log10(umi_counts.clamp(min=eps))  # [cells]

    # Use pd.qcut to get bin indices and edges
    bin_indices_np, bin_edges_np = pd.qcut(
        log_umi.cpu().numpy(), q=n_bins, labels=False, retbins=True, duplicates='drop'
    )

    dtype = gene_detects.dtype
    bin_edges = torch.tensor(bin_edges_np, device=device, dtype=dtype)
    bin_indices = torch.tensor(bin_indices_np, device=device, dtype=dtype)

    # Compute bin centers from bin_edges
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # [n_bins]

    # Compute bin means per gene
    detect_sum = torch.zeros(gene_detects.size(0), len(bin_centers), device=device, dtype=dtype)
    bin_counts = torch.zeros(len(bin_centers), device=device, dtype=dtype)

    for i in range(len(bin_centers)):
        mask = bin_indices == i
        bin_counts[i] = mask.sum()
        if mask.any():
            detect_sum[:, i] = gene_detects[:, mask].sum(dim=1)

    # Laplace smoothing
    bin_detect_rate = (detect_sum + 1) / (bin_counts[None, :] + 2)

    # Fit logit model: y = b0 + b1 * bin_center (per gene)
    logit_y = torch.log(bin_detect_rate / (1 - bin_detect_rate + eps) + eps)  # [genes, bins]
    X = torch.stack([torch.ones_like(bin_centers), bin_centers])  # [2, bins]
    X = X.T  # [bins, 2]

    # Solve (X^T X)^-1 X^T y for each gene
    XT_X = X.T @ X  # [2, 2]
    XT_X_inv = torch.inverse(XT_X)  # [2, 2]
    XT = X.T  # [2, bins]
    b = XT_X_inv @ (XT @ logit_y.T)  # [2, genes]

    b0, b1 = b[0], b[1]  # [genes]
    logit_pred = b0[:, None] + b1[:, None] * log_umi[None, :]  # [genes, cells]
    detect_p = torch.sigmoid(logit_pred)  # [genes, cells]

    mu = detect_p
    var = detect_p * (1 - detect_p)
    x2 = detect_p
    
    return mu, var, x2


def ct_bernoulli_model_linear(gene_detects, umi_counts, cell_types):
    
    mu_ct = np.zeros(len(cell_types))
    var_ct = np.zeros(len(cell_types))
    x2_ct = np.zeros(len(cell_types))
        
    for cell_type in np.unique(cell_types):
        
        gene_detects_ct = gene_detects[cell_types == cell_type]
        umi_counts_ct = umi_counts[cell_types == cell_type]

        # We modify the 0 UMI counts to 1e-10 to remove the NaN values from the qcut output.
        umi_counts_ct[umi_counts_ct == 0] = 1e-10
        
        umi_count_bins, bins = pd.qcut(
            np.log10(umi_counts_ct), N_BIN_TARGET, labels=False, retbins=True,
            duplicates='drop'
        )
        bin_centers = np.array(
            [bins[i] / 2 + bins[i + 1] / 2 for i in range(len(bins) - 1)]
        )

        N_BIN = len(bin_centers)

        bin_detects = bin_gene_detection(gene_detects_ct, umi_count_bins, N_BIN)

        lbin_detects = logit(bin_detects)

        X = np.ones((N_BIN, 2))
        X[:, 1] = bin_centers
        Y = lbin_detects

        b = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
        detect_p = ilogit(b[0] + b[1] * np.log10(umi_counts))

        mu = detect_p
        var = detect_p * (1 - detect_p)
        x2 = detect_p
        
        mu_ct[cell_types == cell_type] = mu
        var_ct[cell_types == cell_type] = var
        x2_ct[cell_types == cell_type] = x2

    return mu_ct, var_ct, x2_ct


bernoulli_model_torch = bernoulli_model_linear_torch
bernoulli_model = bernoulli_model_linear


@njit
def logit(p):
    return np.log(p / (1 - p))


@njit
def ilogit(q):
    return np.exp(q) / (1 + np.exp(q))


@njit
def bin_gene_detection(gene_detects, umi_count_bins, N_BIN):
    bin_detects = np.zeros(N_BIN)
    bin_totals = np.zeros(N_BIN)

    for i in range(len(gene_detects)):
        x = gene_detects[i]
        bin_i = umi_count_bins[i]
        bin_detects[bin_i] += x
        bin_totals[bin_i] += 1

    # Need to account for 0% detects
    #    Add 1 to numerator and denominator
    # Need to account for 100% detects
    #    Add 1 to denominator

    return (bin_detects+1) / (bin_totals+2)


def normal_model(gene_counts, umi_counts):

    """
    Simplest Model - just assumes expression data is normal
    UMI counts are regressed out
    """

    X = np.vstack((np.ones(len(umi_counts)), umi_counts)).T
    y = gene_counts.reshape((-1, 1))

    if umi_counts.var() == 0:

        mu = gene_counts.mean()
        var = gene_counts.var()
        mu = np.repeat(mu, len(umi_counts))
        var = np.repeat(var, len(umi_counts))
        x2 = mu**2 + var

        return mu, var, x2

    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)

    mu = X.dot(B)

    var = (y - mu).var()
    var = np.repeat(var, len(umi_counts))

    mu = mu.ravel()

    x2 = mu**2 + var

    return mu, var, x2


def normal_model_torch(counts: torch.Tensor, umi_counts: torch.Tensor, eps: float = 1e-10):
    """
    Vectorized Normal model in PyTorch for a batch of genes.

    Args:
        counts: Tensor of shape [genes, cells], gene expression values.
        umi_counts: Tensor of shape [cells], total UMI per cell.
        eps: Small constant to avoid instability in matrix inversion.

    Returns:
        mu: Tensor of shape [genes, cells], predicted mean per gene per cell.
        var: Tensor of shape [genes, cells], constant variance per gene.
        x2: Tensor of shape [genes, cells], mu^2 + var.
    """
    device = counts.device
    genes, cells = counts.shape

    # Design matrix X: [cells, 2]
    ones = torch.ones_like(umi_counts).unsqueeze(1)          # [cells, 1]
    umi = umi_counts.unsqueeze(1)                            # [cells, 1]
    X = torch.cat([ones, umi], dim=1)                        # [cells, 2]

    XT = X.T                                                 # [2, cells]
    XT_X = XT @ X                                            # [2, 2]
    try:
        XT_X_inv = torch.inverse(XT_X + eps * torch.eye(2, device=device))  # [2, 2]
    except RuntimeError:
        raise ValueError("Design matrix is singular. Consider regularizing or filtering.")

    # Center y (counts) to [cells, genes] for regression, then transpose
    Y = counts.T                                             # [cells, genes]
    B = XT_X_inv @ (XT @ Y)                                  # [2, genes]

    mu = (X @ B).T                                           # [genes, cells]
    var = ((counts - mu) ** 2).mean(dim=1, keepdim=True)     # [genes, 1]
    var = var.expand_as(mu)                                  # [genes, cells]
    x2 = mu**2 + var                                          # [genes, cells]

    return mu, var, x2


def none_model(gene_counts, umi_counts):

    N = gene_counts.size

    mu = np.zeros(N)
    var = np.ones(N)
    x2 = np.ones(N)

    return mu, var, x2


def none_model_torch(counts: torch.Tensor, umi_counts: torch.Tensor):
    """
    'None' model in PyTorch: returns zero mean and unit variance for all values.

    Args:
        counts: Tensor of shape [genes, cells], ignored here.
        umi_counts: Tensor of shape [cells], ignored here.

    Returns:
        mu: Tensor of zeros [genes, cells]
        var: Tensor of ones [genes, cells]
        x2: Tensor of ones [genes, cells]
    """
    shape = counts.shape
    device = counts.device

    mu = torch.zeros(shape, device=device)
    var = torch.ones(shape, device=device)
    x2 = torch.ones(shape, device=device)

    return mu, var, x2


def apply_model_per_cell_type(
    model_fn: Callable,
    counts: torch.Tensor,
    umi_counts: torch.Tensor,
    cell_types: Union[list, torch.Tensor],
    **kwargs
):
    """
    Applies a model function to each cell type separately.

    Args:
        model_fn: function of form (counts, umi_counts, **kwargs) -> (mu, var, x2)
        counts: [genes, cells] tensor
        umi_counts: [cells] tensor
        cell_types: list or tensor of cell type labels, length = cells
        kwargs: other model-specific arguments

    Returns:
        mu, var, x2: [genes, cells] tensors, concatenated across all cell types
    """
    device = counts.device

    unique_types = cell_types.unique()
    genes, cells = counts.shape

    mu_all = torch.empty((genes, cells), dtype=torch.float64, device=device)
    var_all = torch.empty((genes, cells), dtype=torch.float64, device=device)
    x2_all = torch.empty((genes, cells), dtype=torch.float64, device=device)

    cell_index = np.arange(cells)

    for ct in unique_types:
        
        idx_array = cell_index[cell_types.values == ct]
        idx = torch.tensor(idx_array, device=device)

        counts_ct = counts[:, idx]
        umi_ct = umi_counts[idx]

        mu, var, x2 = model_fn(counts_ct, umi_ct, **kwargs)

        mu_all[:, idx] = mu
        var_all[:, idx] = var
        x2_all[:, idx] = x2

    return mu_all, var_all, x2_all

