import time
from typing import Literal

import numpy as np
import pandas as pd
import pooch
import torch
from anndata import AnnData
from numba import njit
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from scvi.external.harreman._data import harreman_data_hash, harreman_data_url
from scvi.external.harreman.preprocessing.anndata import counts_from_anndata
from scvi.external.harreman.tools.knn import make_weights_non_redundant
from scvi.model._utils import parse_device_args

from . import models

METABOLIC_GENES_DIR = "metabolic_genes"


def load_metabolic_genes(
    species: Literal["mouse"] | Literal["human"] | None = None,
):
    """
    Load the list of metabolic genes for a given species.

    Parameters
    ----------
    species : {"mouse", "human"}, optional (default: "mouse")
        Species used to select the correct metabolic gene list .

    Returns
    -------
    List of metabolic genes.
    """
    fname = f"{species}_metabolic_genes.csv"
    metabolic_genes_path = pooch.retrieve(
        url=harreman_data_url(METABOLIC_GENES_DIR, fname),
        known_hash=harreman_data_hash(METABOLIC_GENES_DIR, fname),
        fname=fname,
        path=pooch.os_cache("scvi_harreman"),
        progressbar=False,
    )

    metabolic_genes = pd.read_csv(metabolic_genes_path, index_col=0)["0"].tolist()

    return metabolic_genes


def compute_local_autocorrelation(
    adata: AnnData,
    layer_key: Literal["use_raw"] | str | None = None,
    database_varm_key: str | None = None,
    model: str | None = None,
    genes: list | None = None,
    use_metabolic_genes: bool = False,
    species: Literal["mouse"] | Literal["human"] | None = "mouse",
    umi_counts_obs_key: str | None = None,
    permutation_test: bool = False,
    M: int | None = 1000,
    seed: int | None = 42,
    check_analytic_null: bool = False,
    device: int | str | torch.device = "auto",
    accelerator: str = "auto",
    verbose: bool = False,
):
    """
    Compute gene-level spatial autocorrelation using spatial weights and centered expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData). Requires `obsp["weights"]` for the spatial graph.
    layer_key : str or "use_raw", optional
        Key in `adata.layers` to use for expression data. Use "use_raw" to access `adata.raw`.
    database_varm_key : str, optional
        Key in `adata.varm` used for filtering genes that are part of the transporter or
        ligand-receptor database.
    model : str, optional
        Normalization model to use for centering gene expression. Options include "none",
        "normal", "bernoulli", or "danb".
    genes : list, optional
        List of gene names to include in the analysis. If `None`, all genes are used or
        selected via metabolic/pathway filters.
    use_metabolic_genes : bool, optional (default: False)
        If `True`, restricts analysis to metabolic genes as defined for the selected species.
    species : {"mouse", "human"}, optional (default: "mouse")
        Species used to select the correct metabolic gene list if `use_metabolic_genes=True`.
    umi_counts_obs_key : str, optional
        Key in `adata.obs` with total UMI counts per cell. If `None`, inferred from the
        expression matrix.
    permutation_test : bool, optional (default: False)
        Whether to compute an empirical p-value and null distribution by permuting the data.
    M : int, optional (default: 1000)
        Number of permutations to perform if `permutation_test` is True.
    seed : int, optional (default: 42)
        Random seed for permutation reproducibility.
    check_analytic_null : bool, optional (default: False)
        Whether to evaluate Z-scores under an analytic null distribution using permutation
        Z-scores.
    device : int, str, or torch.device, optional (default: "auto")
        Device to use for computation. If `"auto"` or a device index, this is resolved with
        :func:`scvi.model._utils.parse_device_args`.
    accelerator : str, optional (default: "auto")
        Accelerator type to use when resolving `device`.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        The results are stored in `adata.uns["gene_autocorrelation_results"]` as a DataFrame.
    """
    start = time.time()
    if verbose:
        print("Computing local autocorrelation...")

    if isinstance(device, str) and device != "auto":
        try:
            device = torch.device(device)
        except (RuntimeError, ValueError):
            pass

    if not isinstance(device, torch.device):
        try:
            _, _, device = parse_device_args(
                accelerator=accelerator,
                devices=device,
                return_device="torch",
                validate_single_device=True,
            )
        except (RuntimeError, ValueError):
            device = torch.device("cpu")
    # Verify the device is actually usable; fall back to CPU if not
    if device.type == "cuda":
        try:
            torch.tensor([0.0], device=device)
        except RuntimeError:
            device = torch.device("cpu")

    adata.uns["layer_key"] = layer_key
    adata.uns["model"] = model
    adata.uns["species"] = species

    sample_specific = "sample_key" in adata.uns

    # Gene selection
    if use_metabolic_genes and genes is None:
        genes = pd.Index(load_metabolic_genes(species)).intersection(adata.var_names)
    elif database_varm_key is not None and genes is None:
        source = adata.raw if (layer_key == "use_raw") else adata
        metab_matrix = source.varm[database_varm_key]
        genes = metab_matrix.loc[(metab_matrix != 0).any(axis=1)].index
    elif genes is None:
        genes = adata.raw.var.index if layer_key == "use_raw" else adata.var_names
    else:
        genes = pd.Index(genes)

    # Load counts
    counts = counts_from_anndata(adata[:, genes], layer_key, dense=True)

    # Gene filtering
    if sample_specific:
        sample_key = adata.uns["sample_key"]
        sample_arr = adata.obs[sample_key].to_numpy()
        mask = np.zeros(counts.shape[0], dtype=bool)
        for sample in np.unique(sample_arr):
            sample_idx = np.where(sample_arr == sample)[0]
            mask |= np.all(counts[:, sample_idx] == 0, axis=1)
    else:
        mask = np.all(counts == 0, axis=1)

    counts = counts[~mask]
    genes = genes[~mask]

    # UMI counts
    num_umi = (
        counts.sum(axis=0)
        if umi_counts_obs_key is None
        else np.asarray(adata.obs[umi_counts_obs_key])
    )
    adata.uns["umi_counts"] = num_umi

    # Convert to tensors
    num_umi = torch.tensor(adata.uns["umi_counts"], dtype=torch.float64, device=device)
    counts = torch.tensor(counts, dtype=torch.float64, device=device)

    # Center values
    counts = standardize_counts(adata, counts, model, num_umi, sample_specific)

    adata.var["local_autocorrelation"] = False
    adata.var.loc[genes, "local_autocorrelation"] = True

    # Compute weights
    weights = make_weights_non_redundant(adata.obsp["weights"]).tocoo()
    Wtot2 = torch.tensor((weights.data**2).sum(), device=device)
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

    # Autocorrelation
    WXt = torch.sparse.mm(weights, counts.T)
    G = (counts.T * WXt).sum(dim=0)
    G_max = 0.5 * torch.sum((counts**2) * D[None, :], dim=1)

    # Results
    results = compute_gene_autocorrelation_results(
        counts=counts,
        weights=weights,
        G=G,
        G_max=G_max,
        Wtot2=Wtot2,
        genes=genes,
        D=D,
        M=M,
        permutation_test=permutation_test,
        seed=seed,
        check_analytic_null=check_analytic_null,
        device=device,
    )

    # Save results
    if isinstance(results, tuple):
        results_df, zs_perm, pvals_perm = results
        adata.uns["analytic_null_ac_zs_perm"] = zs_perm
        adata.uns["analytic_null_ac_pvals_perm"] = pvals_perm
        if verbose:
            print(
                "Analytic null results are stored in adata.uns with the following keys: "
                "['analytic_null_ac_zs_perm', 'analytic_null_ac_pvals_perm']"
            )
    else:
        results_df = results

    results_df = results_df.sort_values("Z", ascending=False)
    results_df.index.name = "Gene"
    cols = ["C", "Z", "Z_Pval", "Z_FDR"]
    if "Perm_Pval" in results_df.columns:
        cols += ["Perm_Pval", "Perm_FDR"]
    adata.uns["gene_autocorrelation_results"] = results_df[cols]
    if verbose:
        print(
            "Local autocorrelation results are stored in adata.uns['gene_autocorrelation_results']"
        )

        print("Finished computing local autocorrelation in %.3f seconds" % (time.time() - start))

    return


def compute_gene_autocorrelation_results(
    counts,
    weights,
    G,
    G_max,
    Wtot2,
    genes,
    D,
    M,
    permutation_test,
    seed,
    check_analytic_null,
    device,
):
    """Compute autocorrelation results from spatial statistics and optional permutations."""
    # Compute core stats
    stats = compute_autocor_Z_scores_torch(G, G_max, Wtot2)

    # Build DataFrame
    results = pd.DataFrame({k: v.cpu().numpy() for k, v in stats.items()}, index=genes)

    # Z P-values and FDR
    results["Z_Pval"] = norm.sf(results["Z"])
    results["Z_FDR"] = multipletests(results["Z_Pval"], method="fdr_bh")[1]

    # Permutation test
    if permutation_test:
        n_genes, n_cells = counts.shape
        perm_array = torch.zeros((n_genes, M), dtype=torch.float16, device=device)

        if check_analytic_null:
            ac_zs_perm_array = torch.zeros((n_genes, M), dtype=torch.float16, device=device)
            ac_pvals_perm_array = torch.zeros((n_genes, M), dtype=torch.float16, device=device)

        torch.manual_seed(seed)
        for i in tqdm(range(M), desc="Permutation test"):
            idx = torch.randperm(n_cells, device=device)
            X_perm = counts[:, idx]  # (genes x cells)
            WXt_perm = torch.sparse.mm(weights, X_perm.T)  # (cells x genes)
            G_perm = (X_perm.T * WXt_perm).sum(dim=0)  # (genes,)
            perm_array[:, i] = G_perm.half()

            if check_analytic_null:
                # Compute G_max for permuted data
                G_perm_max = 0.5 * torch.sum(X_perm**2 * D.unsqueeze(0), dim=1)
                stats_perm = compute_autocor_Z_scores_torch(G_perm, G_perm_max, Wtot2)
                ac_zs_perm_array[:, i] = stats_perm["Z"].half()
                ac_pvals_perm_array[:, i] = torch.tensor(
                    norm.sf(stats_perm["Z"].cpu().numpy()), device=device
                ).half()

        # Compute empirical permutation p-values
        G_expanded = G.unsqueeze(1)  # (genes, 1)
        x = torch.sum(perm_array > G_expanded, dim=1)
        pvals = ((x + 1) / (M + 1)).cpu().numpy()
        results["Perm_Pval"] = pvals
        results["Perm_FDR"] = multipletests(pvals, method="fdr_bh")[1]

        # Save optional nulls
        if check_analytic_null:
            return results, ac_zs_perm_array.cpu().numpy(), ac_pvals_perm_array.cpu().numpy()

    return results


@njit
def center_values(vals, mu, var):
    """Center vals by mu and normalize by the square root of var element-wise."""
    out = np.zeros_like(vals)

    for i in range(len(vals)):
        std = var[i] ** 0.5
        if std == 0:
            out[i] = 0
        else:
            out[i] = (vals[i] - mu[i]) / std

    return out


def center_counts_torch(counts, num_umi, model):
    """Center and normalize counts using the specified model.

    Parameters
    ----------
    counts : torch.Tensor of shape [genes, cells]
    num_umi : torch.Tensor of shape [cells]
    model : str
        One of 'bernoulli', 'danb', 'normal', or 'none'.

    Returns
    -------
    torch.Tensor of shape [genes, cells]
        Centered counts.
    """
    # Binarize if using Bernoulli
    if model == "bernoulli":
        counts = (counts > 0).double()
        mu, var, _ = models.bernoulli_model_torch(counts, num_umi)
    elif model == "danb":
        mu, var, _ = models.danb_model_torch(counts, num_umi)
    elif model == "normal":
        mu, var, _ = models.normal_model_torch(counts, num_umi)
    elif model == "none":
        mu, var, _ = models.none_model_torch(counts, num_umi)
    else:
        raise ValueError(f"Unsupported model type: {model}")

    # Avoid division by zero
    std = torch.sqrt(var)
    std[std == 0] = 1.0

    centered = (counts - mu) / std
    centered[centered == 0] = 0  # Optional: to match old behavior

    return centered


def compute_autocor_Z_scores_torch(G, G_max, Wtot2):
    """Compute autocorrelation Z-scores from G, G_max, and Wtot2.

    Returns a dict with tensors: G, G_max, EG, stdG, Z, C.
    """
    EG = 0.0
    stdG = (Wtot2 - EG**2) ** 0.5

    Z = (G - EG) / stdG  # (genes,)
    C = (G - EG) / G_max  # (genes,)

    EG_tensor = torch.full_like(G, EG)
    stdG_tensor = torch.full_like(G, stdG)

    return {
        "G": G,
        "G_max": G_max,
        "EG": EG_tensor,
        "stdG": stdG_tensor,
        "Z": Z,
        "C": C,
    }


def standardize_counts(adata, counts, model, num_umi, sample_specific):
    """Center and normalize counts using the specified model, optionally per sample."""
    if sample_specific:
        sample_key = adata.uns["sample_key"]
        for sample in adata.obs[sample_key].unique():
            subset = np.where(adata.obs[sample_key] == sample)[0]
            counts[:, subset] = center_counts_torch(counts[:, subset], num_umi[subset], model)
    else:
        counts = center_counts_torch(counts, num_umi, model)

    return counts
