import ast
import itertools
import time
import warnings
from collections import defaultdict
from functools import partial
from typing import Literal

import numpy as np
import pandas as pd
import sparse
import torch
from anndata import AnnData
from numba import jit, njit
from scipy.stats import mannwhitneyu, norm, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from scvi.external.harreman.hotspot import models
from scvi.external.harreman.preprocessing.anndata import counts_from_anndata
from scvi.external.harreman.tools.knn import make_weights_non_redundant


def _lazy_import_hotspot():
    """Resolve circular imports lazily."""
    global compute_local_autocorrelation, standardize_counts
    from scvi.external.harreman.hotspot.local_autocorrelation import (
        compute_local_autocorrelation as _cla,
    )
    from scvi.external.harreman.hotspot.local_autocorrelation import standardize_counts as _sc

    compute_local_autocorrelation = _cla
    standardize_counts = _sc


compute_local_autocorrelation = None
standardize_counts = None


def apply_gene_filtering(
    adata: AnnData,
    layer_key: Literal["use_raw"] | str | None = None,
    cell_type_key: str | None = None,
    model: str | None = None,
    feature_elimination: bool | None = False,
    threshold: float | None = 0.2,
    autocorrelation_filt: bool | None = False,
    expression_filt: bool | None = False,
    de_filt: bool | None = False,
    umi_counts_obs_key: str | None = None,
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: bool | None = False,
):
    """
    Applies multi-step gene filtering to an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData).
    layer_key : str, optional
        Key to use from `adata.layers` or `"use_raw"` to use `adata.raw.X`.
    cell_type_key : str, optional
        Key in `adata.obs` containing cell type annotations.
    model : str, optional
        Model name for autocorrelation computation.
    feature_elimination : bool, optional (default: False)
        If True, filters genes based on sparsity across all cells.
    threshold : float, optional (default: 0.2)
        Minimum fraction of cells in which the gene must be expressed.
    autocorrelation_filt : str, optional (default: False)
        If True, filters genes based on spatial autocorrelation significance.
    expression_filt : str, optional (default: False)
        If True, filters genes based on expression in each cell type.
    de_filt : str, optional (default: False)
        If True, filters genes based on differential expression between each cell type and the rest.
    umi_counts_obs_key : str, optional
        Key in `adata.obs` with total UMI counts per cell. If `None`, inferred from the expression matrix.
    device : torch.device, optional
        Device to use for computation (e.g., CUDA or CPU). Defaults to GPU if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
    """
    start = time.time()
    if verbose:
        print("Applying gene filtering...")

    adata.uns["autocorrelation_filt"] = autocorrelation_filt
    adata.uns["expression_filt"] = expression_filt
    adata.uns["de_filt"] = de_filt

    db_key = adata.uns["database_varm_key"]

    if feature_elimination:
        perform_feature_elimination(adata, layer_key, db_key, threshold)

    _lazy_import_hotspot()
    if autocorrelation_filt:
        compute_local_autocorrelation(
            adata=adata,
            layer_key=layer_key,
            database_varm_key=db_key,
            model=model,
            umi_counts_obs_key=umi_counts_obs_key,
            device=device,
            verbose=verbose,
        )

    if expression_filt or de_filt:
        if cell_type_key is None:
            cell_type_key = adata.uns.get("cell_type_key")
            if cell_type_key is None:
                raise ValueError('The "cell_type_key" argument needs to be provided.')

        filtered_genes, filtered_genes_ct = filter_genes(
            adata, layer_key, db_key, cell_type_key, expression_filt, de_filt, autocorrelation_filt
        )
        adata.uns["filtered_genes"] = filtered_genes
        adata.uns["filtered_genes_ct"] = filtered_genes_ct

    if verbose:
        print("Finished applying gene filtering in %.3f seconds" % (time.time() - start))

    return


def perform_feature_elimination(adata, layer_key, database_varm_key, threshold):
    """
    Filters out genes that are too sparse across all cells.

    Parameters
    ----------
    adata
        Annotated data object (AnnData).
    layer_key
        Which layer to use (or "use_raw").
    database_varm_key
        Key in `adata.varm` pointing to relevant features to filter.
    threshold
        Minimum fraction of cells in which the gene must be expressed.
    """
    use_raw = layer_key == "use_raw"

    metab_matrix = adata.raw.varm[database_varm_key] if use_raw else adata.varm[database_varm_key]
    genes = metab_matrix.loc[(metab_matrix != 0).any(axis=1)].index

    counts = counts_from_anndata(adata[:, genes], layer_key, dense=True)

    valid_genes = genes[filter_expr_matrix(counts, threshold=threshold)]

    adata.varm[database_varm_key][~adata.var_names.isin(valid_genes)] = 0

    return


def filter_genes(
    adata,
    layer_key,
    database_varm_key,
    cell_type_key,
    expression_filt,
    de_filt,
    autocorrelation_filt,
):
    """
    Applies expression and/or DE filtering per cell type.

    Parameters
    ----------
    adata
        Annotated data object (AnnData).
    layer_key
        Which layer to use (or "use_raw").
    database_varm_key
        Key in `adata.varm` pointing to gene database features.
    cell_type_key
        Key in `adata.obs` with cell type labels.
    expression_filt
        Whether to filter based on expression sparsity in each cell type.
    de_filt
        Whether to filter based on differential expression.
    autocorrelation_filt
        Whether to restrict to spatially autocorrelated genes.

    Returns
    -------
    filtered_genes
        List of genes retained across any cell type.
    filtered_genes_ct
        Dict mapping cell types to their filtered genes.
    """
    if autocorrelation_filt:
        autocor_results = adata.uns["gene_autocorrelation_results"]
        sig_genes = autocor_results.query("Z_FDR < 0.05").index
        if len(sig_genes) == 0:
            raise ValueError("There are no significantly autocorrelated genes.")

    else:
        use_raw = layer_key == "use_raw"
        db = adata.raw.varm[database_varm_key] if use_raw else adata.varm[database_varm_key]
        sig_genes = db.loc[(db != 0).any(axis=1)].index

    counts = counts_from_anndata(adata[:, sig_genes], layer_key, dense=True)

    cell_types = adata.obs[cell_type_key].values
    unique_cts = np.unique(cell_types)
    filtered_genes_ct = {}

    # Precompute masks
    masks = {ct: np.where(cell_types == ct)[0] for ct in unique_cts}
    not_masks = {ct: np.where(cell_types != ct)[0] for ct in unique_cts}

    if expression_filt:
        expr_mask = {ct: filter_expr_matrix(counts[:, masks[ct]], 0.2) for ct in unique_cts}
    if de_filt:
        de_stats = {
            ct: de_threshold(counts[:, masks[ct]], counts[:, not_masks[ct]]) for ct in unique_cts
        }

    filtered_genes = set()
    for ct in unique_cts:
        gene_mask = np.ones(len(sig_genes), dtype=bool)

        if expression_filt:
            gene_mask &= expr_mask[ct]

        if de_filt:
            stat, pval, cd = de_stats[ct]
            fdr = multipletests(pval, method="fdr_bh")[1]
            gene_mask &= (fdr < 0.05) & (cd > 0)

        selected = sig_genes[gene_mask]
        filtered_genes_ct[ct] = selected.tolist()
        filtered_genes.update(selected)

    return sorted(filtered_genes), filtered_genes_ct


def filter_expr_matrix(matrix, threshold):

    return (matrix > 0).sum(axis=1) / matrix.shape[1] >= threshold


@njit(parallel=True)
def cohens_d(x, y):

    out = np.empty(x.shape[0])

    for i in range(x.shape[0]):
        nx, ny = len(x[i]), len(y[i])
        vx, vy = np.var(x[i], ddof=1), np.var(y[i], ddof=1)
        pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
        out[i] = (np.mean(x[i]) - np.mean(y[i])) / pooled if pooled > 0 else 0

    return out


def de_threshold(counts_ct, counts_no_ct):

    stat = np.array(
        [
            mannwhitneyu(counts_ct[i], counts_no_ct[i], alternative="greater").statistic
            for i in range(counts_ct.shape[0])
        ]
    )
    pval = np.array(
        [
            mannwhitneyu(counts_ct[i], counts_no_ct[i], alternative="greater").pvalue
            for i in range(counts_ct.shape[0])
        ]
    )
    cd = cohens_d(counts_ct, counts_no_ct)

    return stat, pval, cd


def compute_gene_pairs(
    adata: AnnData,
    layer_key: Literal["use_raw"] | str | None = None,
    cell_type_key: str | None = None,
    cell_type_pairs: list | None = None,
    ct_specific: bool | None = True,
    fix_ct: Literal["all"] | str | None = None,
    verbose: bool | None = False,
):
    """
    Identifies biologically plausible gene pairs involved in ligand-receptor (LR) signaling or
    metabolite transport based on annotated interaction databases and filtered expression data.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData). Must include:
            - `varm["database"]`: DataFrame indicating gene involvement in interactions.
            - `uns["database"]`: 'LR', 'transporter', or 'both'.
            - `uns["ligand"]`, `uns["receptor"]` for LR pairs if applicable.
            - `uns["metabolite_database"]` and/or `uns["LR_database"]` for pair categorization.
            - `obsp["weights"]`: spatial proximity weights.
    layer_key : str or "use_raw", optional
        Specifies the layer or raw data to use for expression filtering.
    cell_type_key : str, optional
        Key in `adata.obs` indicating cell type annotation.
    cell_type_pairs : list of tuple, optional
        List of tuples specifying cell type pairs to consider. If not provided, all combinations are used.
    ct_specific : bool, optional (default: True)
        If True, restrict gene pair computation to combinations relevant to the given cell type annotations.
    fix_ct : str, optional
        Whether to restrict the cell type pairs to a particular cell type.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        Results are stored in the following keys in `adata.uns`: `lcs`, `lc_zs`, `lc_z_pvals`, and `lc_z_FDR`.
    """
    start = time.time()
    if verbose:
        print("Computing gene pairs...")

    from_value_to_type = {
        "LR": {-1.0: "REC", 1.0: "LIG"},
        "transporter": {-1.0: "IMP", 1.0: "EXP", 2.0: "IMP-EXP"},
    }

    # Setup
    layer_key = layer_key or adata.uns.get("layer_key")
    use_raw = layer_key == "use_raw"
    genes = adata.raw.var.index if use_raw else adata.var_names
    adata.uns["fix_ct"] = fix_ct

    if ct_specific:
        cell_type_key = cell_type_key or adata.uns.get("cell_type_key")
        if cell_type_key is None:
            raise ValueError('Please provide the "cell_type_key" argument.')
        adata.uns["cell_type_key"] = cell_type_key
        cell_types = adata.obs[cell_type_key] if not use_raw else adata.raw.obs[cell_type_key]
        cell_types = cell_types.values.astype(str)

    database = adata.varm["database"]

    heterodimer_info = adata.uns.get("heterodimer_info")
    if heterodimer_info is not None:
        heterodimer_info = heterodimer_info.copy()
        heterodimer_info["Genes"] = heterodimer_info["Genes"].apply(ast.literal_eval)

    # Filters
    autocorrelation_filt = adata.uns.get("autocorrelation_filt", False)
    expression_filt = adata.uns.get("expression_filt", False)
    de_filt = adata.uns.get("de_filt", False)

    if expression_filt or de_filt:
        filtered_genes = adata.uns["filtered_genes"]
        filtered_genes_ct = adata.uns["filtered_genes_ct"]
    elif autocorrelation_filt:
        autocor_results = adata.uns["gene_autocorrelation_results"]
        filtered_genes = autocor_results[autocor_results.Z_FDR < 0.05].index.tolist()
    else:
        filtered_genes = list(genes)

    if not filtered_genes:
        raise ValueError("No genes have passed the filters.")

    filtered_genes_set = set(filtered_genes)
    all_genes_set = set(genes)
    non_sig_genes = list(all_genes_set - filtered_genes_set)

    # Filter out uninformative metabolites
    database.loc[non_sig_genes] = 0
    cols_keep = [
        col
        for col in database.columns
        if (
            (np.unique(database[col]) != 0).sum() > 1
            or database[col][database[col] != 0].unique().tolist() == [2]
        )
    ]
    database = database[cols_keep].copy()
    adata.varm["database"] = database

    if ct_specific and "filtered_genes_ct" not in adata.uns:
        filtered_genes_ct = dict.fromkeys(np.unique(cell_types), filtered_genes)
    else:
        filtered_genes_ct = adata.uns.get("filtered_genes_ct", {})

    weights = adata.obsp["weights"]
    if ct_specific:
        if cell_type_pairs is None:
            cell_type_list = list(filtered_genes_ct)
            if fix_ct:
                cell_type_pairs = list(itertools.product(cell_type_list, repeat=2))
                if fix_ct != "all":
                    cell_type_pairs = [pair for pair in cell_type_pairs if pair[0] == fix_ct]
            else:
                cell_type_pairs = list(itertools.combinations_with_replacement(cell_type_list, 2))
        cell_type_pairs_df = pd.Series(cell_type_pairs)
        valid_mask = cell_type_pairs_df.apply(
            get_interacting_cell_type_pairs, args=(weights, cell_types)
        )
        cell_type_pairs = cell_type_pairs_df[valid_mask].tolist()

    # Setup for result aggregation
    gene_pairs_per_metabolite = {}
    gene_pairs = []
    ct_pairs = []
    gene_pairs_per_ct_pair = {} if ct_specific else None

    if adata.uns["database"] == "both":
        metabolites_set = set(adata.uns["metabolite_database"].Metabolite)
        LR_pairs_set = set(adata.uns["LR_database"].index)

    for metabolite in database.columns:
        metab_genes = database.index[database[metabolite] != 0].tolist()
        if not metab_genes:
            continue

        gene_pairs_per_metabolite[metabolite] = {"gene_pair": [], "gene_type": []}

        if adata.uns["database"] == "both":
            if metabolite in metabolites_set:
                int_type = "transporter"
            elif metabolite in LR_pairs_set:
                int_type = "LR"
            else:
                raise ValueError(
                    'The "metabolite" variable needs to be either a metabolite or a LR pair.'
                )
        else:
            int_type = adata.uns["database"]

        # Build gene pairs
        if int_type == "transporter":
            if (
                heterodimer_info is not None
                and metabolite in heterodimer_info["Metabolite"].values
            ):
                for genes_list in heterodimer_info[heterodimer_info["Metabolite"] == metabolite][
                    "Genes"
                ]:
                    if all(g in metab_genes for g in genes_list):
                        metab_genes = [g for g in metab_genes if g not in genes_list] + [
                            tuple(genes_list)
                        ]

            combos = (
                set(itertools.combinations_with_replacement(metab_genes, 2))
                | set(itertools.permutations(metab_genes, 2))
                if ct_specific
                else set(itertools.combinations_with_replacement(metab_genes, 2))
            )
            all_pairs = [
                (list(x) if isinstance(x, tuple) else x, list(y) if isinstance(y, tuple) else y)
                for x, y in combos
            ]

        else:  # LR
            ligand = adata.uns["ligand"].loc[metabolite].dropna().tolist()
            ligand = ligand[0] if len(ligand) == 1 else ligand
            receptor = adata.uns["receptor"].loc[metabolite].dropna().tolist()
            receptor = receptor[0] if len(receptor) == 1 else receptor
            # if len(ligand) == 0 or len(receptor) == 0:
            if not ligand or not receptor:
                continue
            all_pairs = (
                [(ligand, receptor), (receptor, ligand)] if ct_specific else [(ligand, receptor)]
            )

        # Evaluate gene pairs
        for var1, var2 in all_pairs:

            def extract_val(var, metabolite):
                if isinstance(var, str):
                    return database.at[var, metabolite]
                else:
                    vals = list(
                        {
                            database.at[v, metabolite]
                            for v in var
                            if database.at[v, metabolite] != 0
                        }
                    )
                    return vals[0] if len(vals) == 1 else vals

            val1 = extract_val(var1,metabolite)
            val2 = extract_val(var2,metabolite)

            if not val1 or not val2:
                continue
            if val1 == val2 and val1 in (1.0, -1.0):
                continue

            type1 = from_value_to_type[int_type].get(val1)
            type2 = from_value_to_type[int_type].get(val2)

            gene_pairs_per_metabolite[metabolite]["gene_pair"].append((var1, var2))
            gene_pairs_per_metabolite[metabolite]["gene_type"].append((type1, type2))

            if (var1, var2) not in gene_pairs:
                gene_pairs.append((var1, var2))
                if ct_specific:
                    for ct1, ct2 in cell_type_pairs:
                        in_ct1 = (
                            var1 in filtered_genes_ct[ct1]
                            if isinstance(var1, str)
                            else any(v in filtered_genes_ct[ct1] for v in var1)
                        )
                        in_ct2 = (
                            var2 in filtered_genes_ct[ct2]
                            if isinstance(var2, str)
                            else any(v in filtered_genes_ct[ct2] for v in var2)
                        )
                        if in_ct1 and in_ct2:
                            if (ct1, ct2) not in ct_pairs:
                                ct_pairs.append((ct1, ct2))
                            gene_pairs_per_ct_pair.setdefault((ct1, ct2), []).append((var1, var2))

    # Save results
    adata.uns.setdefault("gene_pairs", gene_pairs)
    if ct_specific:
        adata.uns.setdefault("cell_type_pairs", ct_pairs)
        adata.uns.setdefault("gene_pairs_per_ct_pair", gene_pairs_per_ct_pair)
    adata.uns.setdefault("gene_pairs_per_metabolite", gene_pairs_per_metabolite)

    if verbose:
        print("Finished computing gene pairs in %.3f seconds" % (time.time() - start))

    return


def compute_cell_communication(
    adata: AnnData,
    layer_key_p_test: Literal["use_raw"] | str | None = None,
    layer_key_np_test: Literal["use_raw"] | str | None = None,
    model: str = None,
    center_counts_for_np_test: bool | None = False,
    subset_gene_pairs: str | None = None,
    M: int | None = 1000,
    seed: int | None = 42,
    test: Literal["parametric"] | Literal["non-parametric"] | Literal["both"] | None = "both",
    mean: Literal["algebraic"] | Literal["geometric"] | None = "algebraic",
    check_analytic_null: bool | None = False,
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: bool | None = False,
):
    """
    Computes spatially-informed cell-type-agnostic cell-cell communication (CCC) scores and
    significance across all gene pairs using both parametric and non-parametric statistical tests.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. Required fields include:
            - `uns["gene_pairs"]`: list of gene pairs to evaluate.
            - `uns["gene_pairs_per_metabolite"]`: dictionary mapping metabolites to gene pairs.
            - `obsp["weights"]`: sparse matrix encoding spatial cell-cell proximity.
            - (Optional) `uns["LR_database"]`: interaction metadata for pathway annotation.
            - (Optional) `uns["sample_key"]`: if modeling includes sample-specific factors.
    layer_key_p_test : str or "use_raw", optional
        Data layer to use for the parametric test. If `"use_raw"`, uses `adata.raw`.
    layer_key_np_test : str or "use_raw", optional
        Data layer to use for the non-parametric test. If `"use_raw"`, uses `adata.raw`.
    model : str, optional
        Normalization model to use for centering gene expression. Options include "none", "normal", "bernoulli", or "danb".
    center_counts_for_np_test : bool, optional (default: False)
        Whether to center expression counts using the specified model before non-parametric testing.
    subset_gene_pairs : list, optional
        If provided, restricts the analysis to this subset of gene pairs.
    M : int, optional (default: 1000)
        Number of permutations to use if `permutation_test` is True.
    seed : int, optional (default: 42)
        Random seed for permutation reproducibility.
    test : {'parametric', 'non-parametric', 'both'}, optional (default: 'both')
        Specifies which statistical test(s) to run.
    mean : {'algebraic', 'geometric'}, optional (default: 'algebraic')
        Averaging method for multi-gene interactions.
    check_analytic_null : bool, optional (default: False)
        Whether to evaluate Z-scores under an analytic null distribution using permutation Z-scores.
    device : torch.device, optional
        PyTorch device to run computations on. Defaults to CUDA if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        Results are stored in the following `adata.uns` fields:
            - `uns["ccc_results"]["p"]`: Parametric test results (gene pair and metabolite scores, Z, p-values, FDR).
            - `uns["ccc_results"]["np"]`: Non-parametric test results (communication scores, empirical p-values, FDR).
            - `uns["lc_zs"]`: Symmetric matrix of ligand-receptor Z-scores across genes.
            - `uns["gene_pair_dict"]`: Dictionary mapping metabolites to index positions of gene pairs.
            - `uns["D"]`: Vector of total node degrees per cell (spatial connectivity).
            - `uns["genes"]`: Ordered list of involved genes.
            - `uns["gene_pairs_ind"]`: Index-referenced version of `uns["gene_pairs"]`.
    """
    start = time.time()
    if verbose:
        print("Starting cell-cell communication analysis...")

    adata.uns["ccc_results"] = {}

    if test not in ["both", "parametric", "non-parametric"]:
        raise ValueError(
            'The "test" variable should be one of ["both", "parametric", "non-parametric"].'
        )

    if mean not in ["algebraic", "geometric"]:
        raise ValueError('The "mean" variable should be one of ["algebraic", "geometric"].')

    adata.uns["layer_key_p_test"] = layer_key_p_test
    adata.uns["layer_key_np_test"] = layer_key_np_test
    adata.uns["model"] = model
    adata.uns["center_counts_for_np_test"] = center_counts_for_np_test
    adata.uns["mean"] = mean

    run_cell_communication_analysis(
        adata,
        layer_key_p_test,
        layer_key_np_test,
        model,
        center_counts_for_np_test,
        subset_gene_pairs,
        M,
        seed,
        test,
        mean,
        check_analytic_null,
        device,
        verbose,
    )

    if verbose:
        print("Obtaining communication results...")
    get_cell_communication_results(
        adata,
        adata.uns["genes"],
        layer_key_p_test,
        layer_key_np_test,
        model,
        adata.uns["D"],
        test,
        device,
    )

    if verbose:
        print(
            "Finished computing cell-cell communication analysis in %.3f seconds"
            % (time.time() - start)
        )

    return


def run_cell_communication_analysis(
    adata,
    layer_key_p_test,
    layer_key_np_test,
    model,
    center_counts_for_np_test,
    subset_gene_pairs,
    M,
    seed,
    test,
    mean,
    check_analytic_null,
    device,
    verbose,
):

    use_raw = (layer_key_p_test == "use_raw") & (layer_key_np_test == "use_raw")

    cells = (
        adata.raw.obs.index.values.astype(str) if use_raw else adata.obs_names.values.astype(str)
    )

    sample_specific = "sample_key" in adata.uns

    gene_pairs = adata.uns["gene_pairs"] if subset_gene_pairs is None else subset_gene_pairs
    genes = list(np.unique(list(flatten(adata.uns["gene_pairs"]))))
    adata.uns["genes"] = genes
    adata.uns["cells"] = cells

    # Map gene_pairs to index
    gene_pairs_ind = []
    for pair in gene_pairs:
        idx1 = (
            [genes.index(g) for g in pair[0]]
            if isinstance(pair[0], list)
            else genes.index(pair[0])
        )
        idx2 = (
            [genes.index(g) for g in pair[1]]
            if isinstance(pair[1], list)
            else genes.index(pair[1])
        )
        gene_pairs_ind.append((idx1, idx2))
    adata.uns["gene_pairs_ind"] = gene_pairs_ind

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

    adata.uns["D"] = D.cpu().numpy()

    gene_pairs_per_metabolite = adata.uns["gene_pairs_per_metabolite"]

    metabolite_gene_pair_df = pd.DataFrame.from_dict(
        gene_pairs_per_metabolite, orient="index"
    ).reset_index()
    metabolite_gene_pair_df = metabolite_gene_pair_df.rename(columns={"index": "metabolite"})

    metabolite_gene_pair_df["gene_pair"] = metabolite_gene_pair_df["gene_pair"].apply(
        lambda arr: [(sub_array[0], sub_array[1]) for sub_array in arr]
    )
    metabolite_gene_pair_df["gene_type"] = metabolite_gene_pair_df["gene_type"].apply(
        lambda arr: [(sub_array[0], sub_array[1]) for sub_array in arr]
    )

    metabolite_gene_pair_df = pd.concat(
        [
            metabolite_gene_pair_df["metabolite"],
            metabolite_gene_pair_df.explode("gene_pair")["gene_pair"],
            metabolite_gene_pair_df.explode("gene_type")["gene_type"],
        ],
        axis=1,
    )
    metabolite_gene_pair_df = metabolite_gene_pair_df.reset_index(drop=True)

    if "LR_database" in adata.uns.keys():
        LR_database = adata.uns["LR_database"]
        df_merged = pd.merge(
            metabolite_gene_pair_df,
            LR_database,
            left_on="metabolite",
            right_on="interaction_name",
            how="left",
        )
        LR_df = df_merged.dropna(subset=["pathway_name"])
        metabolite_gene_pair_df["metabolite"][
            metabolite_gene_pair_df.metabolite.isin(LR_df.metabolite)
        ] = LR_df["pathway_name"]

    gene_pair_dict = {}
    for metabolite, group in metabolite_gene_pair_df.groupby("metabolite"):
        idxs = (
            group["gene_pair"]
            .apply(lambda gp: gene_pairs.index(gp) if gp in gene_pairs else None)
            .dropna()
            .tolist()
        )
        idxs = [int(ind) for ind in idxs if ind is not None]
        if idxs:
            gene_pair_dict[metabolite] = idxs

    adata.uns["gene_pair_dict"] = gene_pair_dict

    if test in ["parametric", "both"]:
        if verbose:
            print("Running the parametric test...")

        adata.uns["ccc_results"]["p"] = {"gp": {}, "m": {}}

        Wtot2 = torch.tensor((weights.data**2).sum(), device=device)

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_ind:
            if isinstance(idx1, list):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, list):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        # Standardize counts
        _lazy_import_hotspot()
        counts_1 = standardize_counts(adata, counts_1, model, num_umi, sample_specific)
        _lazy_import_hotspot()
        counts_2 = standardize_counts(adata, counts_2, model, num_umi, sample_specific)

        # Compute CCC scores
        WX2t = torch.sparse.mm(weights, counts_2.T)
        WtX2t = torch.sparse.mm(weights.transpose(0, 1), counts_2.T)
        cs_gp = (counts_1.T * WX2t).sum(0) + (counts_1.T * WtX2t).sum(0)
        same_gene_mask = torch.tensor([g1 == g2 for g1, g2 in gene_pairs], device=device)
        cs_gp[same_gene_mask] = cs_gp[same_gene_mask] / 2
        adata.uns["ccc_results"]["p"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()

        # Compute metabolite-level scores
        cs_m = compute_metabolite_cs(cs_gp, gene_pair_dict, interacting_cell_scores=False)
        adata.uns["ccc_results"]["p"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        # Compute second moments
        WX1t = torch.sparse.mm(weights, counts_1.T)
        WtX1t = torch.sparse.mm(weights.transpose(0, 1), counts_1.T)
        eg2_a = (WX1t + WtX1t).pow(2).sum(dim=0)
        eg2_b = (WX2t + WtX2t).pow(2).sum(dim=0)
        eg2s_gp = (eg2_a, eg2_b)

        # Z-score computation
        Z_gp, Z_m = compute_p_results(cs_gp, cs_m, gene_pairs_ind, Wtot2, eg2s_gp, gene_pair_dict)
        # Convert tensors to numpy for statsmodels and pandas
        Z_gp_np = Z_gp.detach().cpu().numpy()
        Z_m_np = Z_m.detach().cpu().numpy()
        # Compute p-values and FDRs
        Z_pvals_gp = norm.sf(Z_gp_np)
        Z_pvals_m = norm.sf(Z_m_np)
        FDR_gp = multipletests(Z_pvals_gp, method="fdr_bh")[1]
        FDR_m = multipletests(Z_pvals_m, method="fdr_bh")[1]

        # Store in AnnData
        adata.uns["ccc_results"]["p"]["gp"]["Z"] = Z_gp_np
        adata.uns["ccc_results"]["p"]["gp"]["Z_pval"] = Z_pvals_gp
        adata.uns["ccc_results"]["p"]["gp"]["Z_FDR"] = FDR_gp
        adata.uns["ccc_results"]["p"]["m"]["Z"] = Z_m_np
        adata.uns["ccc_results"]["p"]["m"]["Z_pval"] = Z_pvals_m
        adata.uns["ccc_results"]["p"]["m"]["Z_FDR"] = FDR_m

        # Symmetric LC Z-score matrix
        genes_ = [
            tuple(i) if isinstance(i, list) else i
            for i in pd.Series([g for pair in gene_pairs for g in pair]).drop_duplicates()
        ]
        gene_pairs_ = [
            (tuple(a) if isinstance(a, list) else a, tuple(b) if isinstance(b, list) else b)
            for a, b in gene_pairs
        ]
        lc_zs = pd.DataFrame(np.zeros((len(genes_), len(genes_))), index=genes_, columns=genes_)
        for i, (g1, g2) in enumerate(gene_pairs_):
            lc_zs.iloc[genes_.index(g1), genes_.index(g2)] = Z_gp_np[i]
        # Force diagonal to 0 and symmetrize
        np.fill_diagonal(lc_zs.values, 0)
        adata.uns["lc_zs"] = (lc_zs + lc_zs.T) / 2

        if verbose:
            print("Parametric test finished.")

    if test in ["non-parametric", "both"]:
        if verbose:
            print("Running the non-parametric test...")

        adata.uns["ccc_results"]["np"] = {"gp": {}, "m": {}}

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_ind:
            if isinstance(idx1, list):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, list):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        if center_counts_for_np_test:
            num_umi = counts.sum(dim=0)
            _lazy_import_hotspot()
            counts_1 = standardize_counts(adata, counts_1, model, num_umi, sample_specific)
            _lazy_import_hotspot()
            counts_2 = standardize_counts(adata, counts_2, model, num_umi, sample_specific)

        n_cells = counts_1.shape[1]
        same_gene_mask = torch.tensor([g1 == g2 for g1, g2 in gene_pairs], device=device)

        if center_counts_for_np_test and test == "both":
            adata.uns["ccc_results"]["np"]["gp"]["cs"] = np.array(
                adata.uns["ccc_results"]["p"]["gp"]["cs"]
            )
            adata.uns["ccc_results"]["np"]["m"]["cs"] = np.array(
                adata.uns["ccc_results"]["p"]["m"]["cs"]
            )
        else:
            WX2t = torch.sparse.mm(weights, counts_2.T)
            WtX2t = torch.sparse.mm(weights.transpose(0, 1), counts_2.T)
            cs_gp = (counts_1.T * WX2t).sum(0) + (counts_1.T * WtX2t).sum(0)
            cs_gp[same_gene_mask] = cs_gp[same_gene_mask] / 2
            adata.uns["ccc_results"]["np"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()
            cs_m = compute_metabolite_cs(cs_gp, gene_pair_dict, interacting_cell_scores=False)
            adata.uns["ccc_results"]["np"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        perm_cs_gp_a = torch.zeros((counts_1.shape[0], M), dtype=torch.float64, device=device)
        perm_cs_gp_b = torch.zeros_like(perm_cs_gp_a)
        perm_cs_m_a = torch.zeros((len(gene_pair_dict), M), dtype=torch.float64, device=device)
        perm_cs_m_b = torch.zeros_like(perm_cs_m_a)

        if check_analytic_null:
            gp_zs_perm_array = torch.zeros_like(perm_cs_gp_a)
            gp_pvals_perm_array = torch.zeros_like(perm_cs_gp_a)
            m_zs_perm_array = torch.zeros_like(perm_cs_m_a)
            m_pvals_perm_array = torch.zeros_like(perm_cs_m_a)

        torch.manual_seed(seed)
        for i in tqdm(range(M), desc="Permutation test"):
            idx = torch.randperm(n_cells, device=device)

            c1_perm_a = counts_1.clone()
            c2_perm_a = counts_2[:, idx]
            c1_perm_a[same_gene_mask] = counts_1[same_gene_mask, :][:, idx]

            WX2t_a = torch.sparse.mm(weights, c2_perm_a.T)
            WtX2t_a = torch.sparse.mm(weights.transpose(0, 1), c2_perm_a.T)
            cs_a = (c1_perm_a.T * WX2t_a).sum(0) + (c1_perm_a.T * WtX2t_a).sum(0)
            cs_a[same_gene_mask] = cs_a[same_gene_mask] / 2
            perm_cs_gp_a[:, i] = cs_a

            cs_m_a = compute_metabolite_cs(cs_a, gene_pair_dict, interacting_cell_scores=False)
            perm_cs_m_a[:, i] = cs_m_a

            c2_perm_b = counts_2.clone()
            c1_perm_b = counts_1[:, idx]
            c2_perm_b[same_gene_mask] = counts_2[same_gene_mask, :][:, idx]

            WX2t_b = torch.sparse.mm(weights, c2_perm_b.T)
            WtX2t_b = torch.sparse.mm(weights.transpose(0, 1), c2_perm_b.T)
            cs_b = (c1_perm_b.T * WX2t_b).sum(0) + (c1_perm_b.T * WtX2t_b).sum(0)
            cs_b[same_gene_mask] = cs_b[same_gene_mask] / 2
            perm_cs_gp_b[:, i] = cs_b

            cs_m_b = compute_metabolite_cs(cs_b, gene_pair_dict, interacting_cell_scores=False)
            perm_cs_m_b[:, i] = cs_m_b

            if check_analytic_null:
                Z_gp_perm, Z_m_perm = compute_p_results(
                    (cs_a, cs_b), (cs_m_a, cs_m_b), gene_pairs_ind, Wtot2, eg2s_gp, gene_pair_dict
                )
                gp_zs_perm_array[:, i] = Z_gp_perm
                gp_pvals_perm_array[:, i] = torch.tensor(
                    norm.sf(Z_gp_perm.cpu().numpy()), device=device
                )
                m_zs_perm_array[:, i] = Z_m_perm
                m_pvals_perm_array[:, i] = torch.tensor(
                    norm.sf(Z_m_perm.cpu().numpy()), device=device
                )

        adata.uns["ccc_results"]["np"]["gp"]["perm_cs_a"] = perm_cs_gp_a.detach().cpu().numpy()
        adata.uns["ccc_results"]["np"]["gp"]["perm_cs_b"] = perm_cs_gp_b.detach().cpu().numpy()
        adata.uns["ccc_results"]["np"]["m"]["perm_cs_a"] = perm_cs_m_a.detach().cpu().numpy()
        adata.uns["ccc_results"]["np"]["m"]["perm_cs_b"] = perm_cs_m_b.detach().cpu().numpy()

        x_gp_a = (perm_cs_gp_a > cs_gp[:, None]).sum(dim=1)
        x_gp_b = (perm_cs_gp_b > cs_gp[:, None]).sum(dim=1)
        x_m_a = (perm_cs_m_a > cs_m[:, None]).sum(dim=1)
        x_m_b = (perm_cs_m_b > cs_m[:, None]).sum(dim=1)

        pvals_gp_a = (x_gp_a + 1).float() / (M + 1)
        pvals_gp_b = (x_gp_b + 1).float() / (M + 1)
        pvals_m_a = (x_m_a + 1).float() / (M + 1)
        pvals_m_b = (x_m_b + 1).float() / (M + 1)

        pvals_gp = torch.where(pvals_gp_a > pvals_gp_b, pvals_gp_a, pvals_gp_b)
        pvals_m = torch.where(pvals_m_a > pvals_m_b, pvals_m_a, pvals_m_b)

        adata.uns["ccc_results"]["np"]["gp"]["pval"] = pvals_gp.cpu().numpy()
        adata.uns["ccc_results"]["np"]["gp"]["FDR"] = multipletests(
            pvals_gp.cpu().numpy(), method="fdr_bh"
        )[1]
        adata.uns["ccc_results"]["np"]["m"]["pval"] = pvals_m.cpu().numpy()
        adata.uns["ccc_results"]["np"]["m"]["FDR"] = multipletests(
            pvals_m.cpu().numpy(), method="fdr_bh"
        )[1]

        if check_analytic_null:
            adata.uns["ccc_results"]["np"]["analytic_null"] = {
                "gp_zs_perm": gp_zs_perm_array.detach().cpu().numpy(),
                "gp_pvals_perm": gp_pvals_perm_array.detach().cpu().numpy(),
                "m_zs_perm": m_zs_perm_array.detach().cpu().numpy(),
                "m_pvals_perm": m_pvals_perm_array.detach().cpu().numpy(),
            }

    if verbose:
        print("Non-parametric test finished.")

    return


def compute_ct_cell_communication(
    adata: AnnData,
    layer_key_p_test: Literal["use_raw"] | str | None = None,
    layer_key_np_test: Literal["use_raw"] | str | None = None,
    model: str = None,
    cell_type_key: str | None = None,
    center_counts_for_np_test: bool | None = False,
    subset_gene_pairs: list | None = None,
    subset_metabolites: list | None = None,
    fix_gp: bool | None = False,
    M: int | None = 1000,
    seed: int | None = 42,
    test: Literal["parametric"] | Literal["non-parametric"] | Literal["both"] | None = "both",
    mean: Literal["algebraic"] | Literal["geometric"] | None = "algebraic",
    check_analytic_null: bool | None = False,
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: bool | None = False,
):
    """
    Computes cell type-aware cell-cell communication (CCC) scores by stratifying communication
    by interacting cell type pairs. Supports parametric and non-parametric statistical inference.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. Required fields include:
            - `uns["gene_pairs"]`: gene pairs involved in communication.
            - `uns["gene_pairs_per_metabolite"]`: maps metabolites to gene pairs.
            - `uns["gene_pairs_per_ct_pair"]`: gene pairs per cell type pair.
            - `obsp["weights"]`: sparse cell-cell proximity matrix.
            - `obs[cell_type_key]`: categorical cell type annotations.
            - `uns["cell_type_pairs"]`: list of interacting cell type pairs.
            - (Optional) `uns["LR_database"]`: for metabolite/pathway annotation.
    layer_key_p_test : str or "use_raw", optional
        Data layer to use for parametric test.
    layer_key_np_test : str or "use_raw", optional
        Data layer to use for non-parametric test.
    model : str, optional
        Normalization model to use for centering gene expression. Options include "none", "normal", "bernoulli", or "danb".
    cell_type_key : str, optional
        Key in `adata.obs` corresponding to cell type annotations. Required if not stored in `uns`.
    center_counts_for_np_test : bool, optional (default: False)
        Whether to center expression counts using the specified model before non-parametric testing.
    subset_gene_pairs : list, optional
        Subset of gene pairs to consider. If None, uses all pairs.
    subset_metabolites : list, optional
        Subset of metabolites to include in the analysis.
    fix_gp : bool, optional (default: False)
        If True, keeps gene pair identity fixed during permutation testing, randomizing cell types only.
    M : int, optional (default: 1000)
        Number of permutations to use if `permutation_test` is True.
    seed : int, optional (default: 42)
        Random seed for permutation reproducibility.
    test : {'parametric', 'non-parametric', 'both'}, optional (default: 'both')
        Specifies which statistical test(s) to run.
    mean : {'algebraic', 'geometric'}, optional (default: 'algebraic')
        Averaging method for multi-gene modules.
    check_analytic_null : bool, optional (default: False)
        Whether to compute Z-scores and p-values under the null distribution for the permutation test.
    device : torch.device, optional
        PyTorch device to run computations on. Defaults to CUDA if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        Results are stored in the following `adata.uns` fields:
            - `ct_ccc_results["p"]`: parametric test results (scores, Z, p-values, FDRs) per gene pair and metabolite per cell type pair.
            - `ct_ccc_results["np"]`: non-parametric test results (communication scores, empirical p-values, FDRs).
            - `gene_pair_dict`: dictionary mapping metabolites to relevant gene pairs.
            - `gene_pairs_ind`, `gene_pairs_ind_per_ct_pair`: index-referenced gene pair representations.
            - `D`: spatial node degree for each cell per cell type pair.
            - `cells`, `genes`: ordered list of cells and genes used in analysis.
            - (optional) `ct_ccc_results["np"]["analytic_null"]`: null distributions from permutation test Z-scores and p-values.
    """
    start = time.time()
    if verbose:
        print("Starting cell type-aware cell-cell communication analysis...")

    adata.uns["ct_ccc_results"] = {}

    if test not in ["both", "parametric", "non-parametric"]:
        raise ValueError(
            'The "test" variable should be one of ["both", "parametric", "non-parametric"].'
        )

    if mean not in ["algebraic", "geometric"]:
        raise ValueError('The "mean" variable should be one of ["algebraic", "geometric"].')

    if "cell_type_key" in adata.uns and cell_type_key is None:
        cell_type_key = adata.uns["cell_type_key"]
    elif "cell_type_key" not in adata.uns and cell_type_key is None:
        raise ValueError('Please provide the "cell_type_key" argument.')

    adata.uns["layer_key_p_test"] = layer_key_p_test
    adata.uns["layer_key_np_test"] = layer_key_np_test
    adata.uns["model"] = model
    adata.uns["cell_type_key"] = cell_type_key
    adata.uns["center_counts_for_np_test"] = center_counts_for_np_test
    adata.uns["mean"] = mean

    run_ct_cell_communication_analysis(
        adata,
        layer_key_p_test,
        layer_key_np_test,
        model,
        cell_type_key,
        center_counts_for_np_test,
        subset_gene_pairs,
        subset_metabolites,
        fix_gp,
        M,
        seed,
        test,
        mean,
        check_analytic_null,
        device,
        verbose,
    )

    if verbose:
        print("Obtaining cell type-aware communication results...")
    get_ct_cell_communication_results(
        adata,
        adata.uns["genes"],
        adata.uns["cells"],
        layer_key_p_test,
        layer_key_np_test,
        model,
        adata.obs[cell_type_key],
        adata.uns["cell_type_pairs"],
        adata.uns["D"],
        test,
        device,
    )

    if verbose:
        print(
            "Finished computing cell type-aware cell-cell communication analysis in %.3f seconds"
            % (time.time() - start)
        )

    return


def run_ct_cell_communication_analysis(
    adata,
    layer_key_p_test,
    layer_key_np_test,
    model,
    cell_type_key,
    center_counts_for_np_test,
    subset_gene_pairs,
    subset_metabolites,
    fix_gp,
    M,
    seed,
    test,
    mean,
    check_analytic_null,
    device,
    verbose,
):

    use_raw = (layer_key_p_test == "use_raw") & (layer_key_np_test == "use_raw")
    obs = adata.raw.obs if use_raw else adata.obs
    cells = (
        adata.raw.obs.index.values.astype(str) if use_raw else adata.obs_names.values.astype(str)
    )

    sample_specific = "sample_key" in adata.uns

    fix_ct = True if adata.uns["fix_ct"] else False

    gene_pairs = adata.uns["gene_pairs"] if subset_gene_pairs is None else subset_gene_pairs
    genes = list(np.unique(list(flatten(adata.uns["gene_pairs"]))))
    adata.uns["genes"] = genes

    cell_types = obs[cell_type_key]
    cell_type_pairs = adata.uns.get("cell_type_pairs")
    gene_pairs_per_ct_pair = adata.uns.get("gene_pairs_per_ct_pair", {})

    weights = adata.obsp["weights"]

    used_ct_pairs = list({ct for cell_type_pair in cell_type_pairs for ct in cell_type_pair})
    all_cell_types = set(cell_types.unique())
    used_ct_pairs_set = set(used_ct_pairs)
    if used_ct_pairs_set < all_cell_types:
        keep_mask = cell_types[cells].isin(used_ct_pairs).values
        keep_indices = np.where(keep_mask)[0]
        weights = weights[keep_indices][:, keep_indices]
        cells = cells[keep_indices]
        cell_types = cell_types.loc[cells]

    adata.uns["cells"] = cells

    weights_ct_pairs = create_weights_ct_pairs(
        weights.tocoo(), cell_types, cell_type_pairs, device
    )

    row_degrees = torch.sparse.sum(weights_ct_pairs, dim=2).to_dense()
    col_degrees = torch.sparse.sum(weights_ct_pairs, dim=1).to_dense()
    D = row_degrees + col_degrees
    if used_ct_pairs_set < all_cell_types:
        D_full = torch.zeros(
            (len(cell_type_pairs), adata.shape[0]),
            device=weights_ct_pairs.device,
            dtype=weights_ct_pairs.dtype,
        )
        D_full[:, keep_indices] = D
        adata.uns["D"] = D_full.cpu().numpy()
    else:
        adata.uns["D"] = D.cpu().numpy()

    # Map gene_pairs to index
    gene_pairs_ind = []
    for pair in gene_pairs:
        idx1 = (
            [genes.index(g) for g in pair[0]]
            if isinstance(pair[0], list)
            else genes.index(pair[0])
        )
        idx2 = (
            [genes.index(g) for g in pair[1]]
            if isinstance(pair[1], list)
            else genes.index(pair[1])
        )
        gene_pairs_ind.append((idx1, idx2))
    adata.uns["gene_pairs_ind"] = gene_pairs_ind

    # Cell-type pair-specific indices
    gene_pairs_ind_per_ct_pair = defaultdict(list)
    gene_pairs_per_ct_pair_ind = defaultdict(list)
    for ct_pair, gpairs in gene_pairs_per_ct_pair.items():
        for pair in gpairs:
            if pair not in gene_pairs:
                continue
            idx = gene_pairs.index(pair)
            gene_pairs_ind_per_ct_pair[ct_pair].append(gene_pairs_ind[idx])
            gene_pairs_per_ct_pair_ind[ct_pair].append(idx)

    adata.uns["gene_pairs_ind_per_ct_pair"] = dict(gene_pairs_ind_per_ct_pair)
    adata.uns["gene_pairs_per_ct_pair_ind"] = dict(gene_pairs_per_ct_pair_ind)

    def make_hashable(pair):
        return tuple(tuple(x) if isinstance(x, list) else x for x in pair)

    gene_pairs_ind_set = {make_hashable(pair) for pair in gene_pairs_ind}
    ct_specific_gene_pairs = [
        i
        for i, pairs in enumerate(gene_pairs_ind_per_ct_pair.values())
        if {make_hashable(pair) for pair in pairs} < gene_pairs_ind_set
    ]

    # Metabolite-gene pair preparation
    gp_metab = adata.uns["gene_pairs_per_metabolite"]
    metabolite_gene_pair_df = (
        pd.DataFrame.from_dict(gp_metab, orient="index")
        .rename_axis("metabolite")
        .explode(["gene_pair", "gene_type"])
        .reset_index()
    )

    if "LR_database" in adata.uns:
        merged = metabolite_gene_pair_df.merge(
            adata.uns["LR_database"], left_on="metabolite", right_on="interaction_name", how="left"
        )
        LR_df = merged.dropna(subset=["pathway_name"])
        metabolite_gene_pair_df.loc[
            metabolite_gene_pair_df.metabolite.isin(LR_df.metabolite), "metabolite"
        ] = LR_df["pathway_name"].values

    if subset_metabolites:
        metabolite_gene_pair_df = metabolite_gene_pair_df[
            metabolite_gene_pair_df.metabolite.isin(subset_metabolites)
        ]

    gene_pair_dict = {}
    for metabolite, group in metabolite_gene_pair_df.groupby("metabolite"):
        idxs = (
            group["gene_pair"]
            .apply(lambda gp: gene_pairs.index(gp) if gp in gene_pairs else None)
            .dropna()
            .tolist()
        )
        idxs = [int(ind) for ind in idxs if ind is not None]
        if idxs:
            gene_pair_dict[metabolite] = idxs
    adata.uns["gene_pair_dict"] = gene_pair_dict

    if test in ["parametric", "both"]:
        if verbose:
            print("Running the parametric test...")

        adata.uns["ct_ccc_results"]["p"] = {"gp": {}, "m": {}}

        weights_sq_data = weights_ct_pairs.values() ** 2
        weights_sq = torch.sparse_coo_tensor(
            weights_ct_pairs.indices(),
            weights_sq_data,
            weights_ct_pairs.shape,
            device=weights_ct_pairs.device,
        )
        Wtot2 = torch.sparse.sum(weights_sq, dim=(1, 2)).to_dense()

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_ind:
            if isinstance(idx1, list):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, list):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        counts_1 = standardize_ct_counts(
            adata, counts_1, model, num_umi, sample_specific, cell_types
        )
        counts_2 = standardize_ct_counts(
            adata, counts_2, model, num_umi, sample_specific, cell_types
        )

        # Compute CCC scores
        cs_gp = torch.zeros((len(cell_type_pairs), counts_1.shape[0]), device=counts_1.device)
        for ct_pair in range(len(cell_type_pairs)):
            W = weights_ct_pairs[ct_pair].coalesce()
            WX2t = torch.sparse.mm(W, counts_2.T)
            cs_gp[ct_pair] = (counts_1.T * WX2t).sum(0)
        adata.uns["ct_ccc_results"]["p"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()

        cs_m = compute_metabolite_cs_ct(
            cs_gp,
            cell_type_key,
            gene_pair_dict,
            gene_pairs_per_ct_pair_ind,
            ct_specific_gene_pairs,
            interacting_cell_scores=False,
        )
        adata.uns["ct_ccc_results"]["p"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        EG2_gp = torch.zeros_like(cs_gp) if fix_ct or fix_gp else Wtot2
        if fix_ct:
            for ct_pair in range(len(cell_type_pairs)):
                W = weights_ct_pairs[ct_pair].coalesce()
                W_sq_data = W.values() ** 2
                W_sq = torch.sparse_coo_tensor(W.indices(), W_sq_data, W.shape, device=W.device)
                X1_sq = counts_1**2
                EG2_gp[ct_pair] = torch.sparse.mm(W_sq, X1_sq.T).sum(0)
        elif fix_gp:
            for ct_pair in range(len(cell_type_pairs)):
                W = weights_ct_pairs[ct_pair].coalesce()
                W_sq_data = W.values() ** 2
                W_sq = torch.sparse_coo_tensor(W.indices(), W_sq_data, W.shape, device=W.device)
                X1_sq = counts_1**2
                X2_sq = counts_2**2
                EG2_gp[ct_pair] = (X1_sq.T * torch.sparse.mm(W_sq, X2_sq.T)).sum(0)

        Z_gp, Z_m = compute_ct_p_results(
            cs_gp,
            cs_m,
            gene_pairs_per_ct_pair_ind,
            ct_specific_gene_pairs,
            EG2_gp,
            cell_type_key,
            gene_pair_dict,
        )

        # Convert tensors to numpy for statsmodels and pandas
        Z_gp_np = Z_gp.detach().cpu().numpy()
        Z_m_np = Z_m.detach().cpu().numpy()
        # Compute p-values and FDRs
        Z_pvals_gp = norm.sf(Z_gp_np)
        Z_pvals_m = norm.sf(Z_m_np)
        FDR_gp = multipletests(Z_pvals_gp.flatten(), method="fdr_bh")[1].reshape(Z_pvals_gp.shape)
        FDR_m = multipletests(Z_pvals_m.flatten(), method="fdr_bh")[1].reshape(Z_pvals_m.shape)

        # Store in AnnData
        adata.uns["ct_ccc_results"]["p"]["gp"]["Z"] = Z_gp_np
        adata.uns["ct_ccc_results"]["p"]["gp"]["Z_pval"] = Z_pvals_gp
        adata.uns["ct_ccc_results"]["p"]["gp"]["Z_FDR"] = FDR_gp
        adata.uns["ct_ccc_results"]["p"]["m"]["Z"] = Z_m_np
        adata.uns["ct_ccc_results"]["p"]["m"]["Z_pval"] = Z_pvals_m
        adata.uns["ct_ccc_results"]["p"]["m"]["Z_FDR"] = FDR_m

    if test in ["non-parametric", "both"]:
        if verbose:
            print("Running the non-parametric test...")

        adata.uns["ct_ccc_results"]["np"] = {"gp": {}, "m": {}}

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_ind:
            if isinstance(idx1, list):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, list):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        if center_counts_for_np_test:
            num_umi = counts.sum(dim=0)
            counts_1 = standardize_ct_counts(
                adata, counts_1, model, num_umi, sample_specific, cell_types
            )
            counts_2 = standardize_ct_counts(
                adata, counts_2, model, num_umi, sample_specific, cell_types
            )

        if center_counts_for_np_test and test == "both":
            adata.uns["ct_ccc_results"]["np"]["gp"]["cs"] = np.array(
                adata.uns["ct_ccc_results"]["p"]["gp"]["cs"]
            )
            adata.uns["ct_ccc_results"]["np"]["m"]["cs"] = np.array(
                adata.uns["ct_ccc_results"]["p"]["m"]["cs"]
            )
        else:
            cs_gp = torch.zeros((len(cell_type_pairs), counts_1.shape[0]), device=counts_1.device)
            for ct_pair in range(len(cell_type_pairs)):
                W = weights_ct_pairs[ct_pair].coalesce()
                WX2t = torch.sparse.mm(W, counts_2.T)
                cs_gp[ct_pair] = (counts_1.T * WX2t).sum(0)
            adata.uns["ct_ccc_results"]["np"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()
            cs_m = compute_metabolite_cs_ct(
                cs_gp,
                cell_type_key,
                gene_pair_dict,
                gene_pairs_per_ct_pair_ind,
                ct_specific_gene_pairs,
                interacting_cell_scores=False,
            )
            adata.uns["ct_ccc_results"]["np"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        perm_cs_gp = torch.zeros(
            (len(cell_type_pairs), counts_1.shape[0], M), dtype=torch.float64, device=device
        )
        perm_cs_m = torch.zeros(
            (len(cell_type_pairs), len(gene_pair_dict), M), dtype=torch.float64, device=device
        )

        if check_analytic_null:
            gp_zs_perm_array = torch.zeros_like(perm_cs_gp)
            gp_pvals_perm_array = torch.zeros_like(perm_cs_gp)
            m_zs_perm_array = torch.zeros_like(perm_cs_m)
            m_pvals_perm_array = torch.zeros_like(perm_cs_m)

        if fix_gp:
            c1_perm = counts_1
            c2_perm = counts_2

        torch.manual_seed(seed)
        for i in tqdm(range(M), desc="Permutation test"):
            if fix_gp:
                indices = torch.randperm(len(cell_types)).numpy()
                shuffled_cell_types = cell_types.iloc[indices].reset_index(drop=True)
                weights_ct_pairs = create_weights_ct_pairs(
                    weights.tocoo(), shuffled_cell_types, cell_type_pairs, device
                )
            else:
                cell_type_labels = torch.tensor(
                    cell_types.astype("category").cat.codes.values, device=counts_1.device
                )
                idx = torch.empty_like(cell_type_labels, dtype=torch.int64)

                for ct in torch.unique(cell_type_labels):
                    ct_mask = cell_type_labels == ct
                    ct_indices = torch.nonzero(ct_mask, as_tuple=True)[0]
                    permuted_indices = ct_indices[torch.randperm(len(ct_indices))]
                    idx[ct_indices] = permuted_indices

                c1_perm = counts_1 if fix_ct else counts_1[:, idx.long()]
                c2_perm = counts_2[:, idx.long()]

            cs_gp = torch.zeros((len(cell_type_pairs), c1_perm.shape[0]), device=c1_perm.device)
            for ct_pair in range(len(cell_type_pairs)):
                W = weights_ct_pairs[ct_pair].coalesce()
                WX2t = torch.sparse.mm(W, c2_perm.T)
                cs_gp[ct_pair] = (c1_perm.T * WX2t).sum(0)
            perm_cs_gp[:, :, i] = cs_gp

            cs_m = compute_metabolite_cs_ct(
                cs_gp,
                cell_type_key,
                gene_pair_dict,
                gene_pairs_per_ct_pair_ind,
                ct_specific_gene_pairs,
                interacting_cell_scores=False,
            )
            perm_cs_m[:, :, i] = cs_m

            if check_analytic_null:
                Z_gp_perm, Z_m_perm = compute_ct_p_results(
                    cs_gp,
                    cs_m,
                    gene_pairs_per_ct_pair_ind,
                    ct_specific_gene_pairs,
                    EG2_gp,
                    cell_type_key,
                    gene_pair_dict,
                )
                gp_zs_perm_array[:, :, i] = Z_gp_perm
                gp_pvals_perm_array[:, :, i] = torch.tensor(
                    norm.sf(Z_gp_perm.cpu().numpy()), device=device
                )
                m_zs_perm_array[:, :, i] = Z_m_perm
                m_pvals_perm_array[:, :, i] = torch.tensor(
                    norm.sf(Z_m_perm.cpu().numpy()), device=device
                )

        adata.uns["ct_ccc_results"]["np"]["gp"]["perm_cs"] = perm_cs_gp.detach().cpu().numpy()
        adata.uns["ct_ccc_results"]["np"]["m"]["perm_cs"] = perm_cs_m.detach().cpu().numpy()

        x_gp = np.sum(
            adata.uns["ct_ccc_results"]["np"]["gp"]["perm_cs"]
            > adata.uns["ct_ccc_results"]["np"]["gp"]["cs"][:, :, np.newaxis],
            axis=2,
        )
        x_m = np.sum(
            adata.uns["ct_ccc_results"]["np"]["m"]["perm_cs"]
            > adata.uns["ct_ccc_results"]["np"]["m"]["cs"][:, :, np.newaxis],
            axis=2,
        )

        pvals_gp = (x_gp + 1) / (M + 1)
        pvals_m = (x_m + 1) / (M + 1)

        adata.uns["ct_ccc_results"]["np"]["gp"]["pval"] = pvals_gp
        adata.uns["ct_ccc_results"]["np"]["gp"]["FDR"] = multipletests(
            pvals_gp.flatten(), method="fdr_bh"
        )[1].reshape(pvals_gp.shape)
        adata.uns["ct_ccc_results"]["np"]["m"]["pval"] = pvals_m
        adata.uns["ct_ccc_results"]["np"]["m"]["FDR"] = multipletests(
            pvals_m.flatten(), method="fdr_bh"
        )[1].reshape(pvals_m.shape)

        if check_analytic_null:
            adata.uns["ct_ccc_results"]["np"]["analytic_null"] = {
                "gp_zs_perm": gp_zs_perm_array.detach().cpu().numpy(),
                "gp_pvals_perm": gp_pvals_perm_array.detach().cpu().numpy(),
                "m_zs_perm": m_zs_perm_array.detach().cpu().numpy(),
                "m_pvals_perm": m_pvals_perm_array.detach().cpu().numpy(),
            }

    adata.uns["cell_types"] = cell_types.tolist() if cell_type_key else None

    if verbose:
        print("Non-parametric test finished.")

    return


def standardize_ct_counts(adata, counts, model, num_umi, sample_specific, cell_types):

    if sample_specific:
        sample_key = adata.uns["sample_key"]
        for sample in adata.obs[sample_key].unique():
            subset = np.where(adata.obs[sample_key] == sample)[0]
            counts[:, subset] = center_ct_counts_torch(
                counts[:, subset], num_umi[subset], model, cell_types[subset]
            )
    else:
        counts = center_ct_counts_torch(counts, num_umi, model, cell_types)

    return counts


def flatten(nested_list):
    for item in nested_list:
        if isinstance(item, list | tuple):
            yield from flatten(item)
        else:
            yield item


def create_weights_ct_pairs(weights, cell_types, cell_type_pairs, device):

    indices = torch.tensor([weights.row, weights.col], dtype=torch.long, device=device)
    values = torch.tensor(weights.data, dtype=torch.float64, device=device)
    shape = weights.shape

    cell_type_cats = cell_types.astype("category")
    cell_type_codes = torch.tensor(
        cell_type_cats.cat.codes.values, dtype=torch.long, device=device
    )
    ct_name_to_code = {name: code for code, name in enumerate(cell_type_cats.cat.categories)}

    row_idx, col_idx = indices
    sender_types = cell_type_codes[row_idx]
    receiver_types = cell_type_codes[col_idx]

    weights_list = []
    coord_list = []

    for i, (ct1, ct2) in enumerate(cell_type_pairs):
        code1 = ct_name_to_code[ct1]
        code2 = ct_name_to_code[ct2]

        pair_mask = (sender_types == code1) & (receiver_types == code2)
        if pair_mask.sum() == 0:
            continue

        pair_values = values[pair_mask]
        pair_coords = torch.stack(
            [
                torch.full((pair_values.shape[0],), i, dtype=torch.long, device=device),
                row_idx[pair_mask],
                col_idx[pair_mask],
            ],
            dim=0,
        )

        weights_list.append(pair_values)
        coord_list.append(pair_coords)

    all_values = torch.cat(weights_list)
    all_coords = torch.cat(coord_list, dim=1)
    weights_ct_pairs = torch.sparse_coo_tensor(
        all_coords, all_values, (len(cell_type_pairs), shape[0], shape[1]), device=device
    )
    weights_ct_pairs = weights_ct_pairs.coalesce()

    return weights_ct_pairs


def select_significant_interactions(
    adata: AnnData,
    ct_aware: bool | None = False,
    test: Literal["parametric"] | Literal["non-parametric"] | None = "parametric",
    use_FDR: bool | None = True,
    threshold: float | None = 0.05,
):
    """
    Select significant gene pairs or metabolite-mediated interactions based on
    FDR/p-value thresholds and (optionally) cell-type–aware tests.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing:
        - ``uns['ccc_results']`` or ``uns['ct_ccc_results']``, each of which includes:
            * ``cell_com_df_gp``: DataFrame with statistics for gene pairs
            * ``cell_com_df_m``:  DataFrame with statistics for metabolites
    ct_aware : bool, default False
        If True, use cell-type–aware CCC results (``uns['ct_ccc_results']``).
        If False, use cell-type-agnostic CCC results (``uns['ccc_results']``).
    test : {"parametric", "non-parametric"}, default "parametric"
        Determines which statistical columns to use:
        - Parametric: ``Z_FDR`` / ``Z_pval``, ``C_p``
        - Non-parametric: ``FDR_np`` / ``pval_np``, ``C_np``
    use_FDR : bool, default True
        If True, threshold significance using FDR values.
        If False, use raw p-values.
    threshold : float, default 0.05
        Significance cutoff applied to the selected statistic (FDR or p-value).
    """
    ccc_key = "ct_ccc_results" if ct_aware else "ccc_results"
    sig_key = "FDR" if use_FDR else "pval"

    if test == "parametric":
        FDR_values_gp = adata.uns[ccc_key]["cell_com_df_gp"][f"Z_{sig_key}"].values
        C_values_gp = adata.uns[ccc_key]["cell_com_df_gp"]["C_p"].values
        FDR_values_m = adata.uns[ccc_key]["cell_com_df_m"][f"Z_{sig_key}"].values
        C_values_m = adata.uns[ccc_key]["cell_com_df_m"]["C_p"].values
    elif test == "non-parametric":
        FDR_values_gp = adata.uns[ccc_key]["cell_com_df_gp"][f"{sig_key}_np"].values
        C_values_gp = adata.uns[ccc_key]["cell_com_df_gp"]["C_np"].values
        FDR_values_m = adata.uns[ccc_key]["cell_com_df_m"][f"{sig_key}_np"].values
        C_values_m = adata.uns[ccc_key]["cell_com_df_m"]["C_np"].values
    else:
        raise ValueError('The "test" variable should be one of ["parametric", "non-parametric"].')

    # Gene pair
    adata.uns[ccc_key]["cell_com_df_gp"]["selected"] = (
        (FDR_values_gp < threshold) & (C_values_gp > 0)
        if test == "non-parametric"
        else (FDR_values_gp < threshold)
    )
    cell_com_df_gp = adata.uns[ccc_key]["cell_com_df_gp"]
    adata.uns[ccc_key]["cell_com_df_gp_sig"] = cell_com_df_gp[
        cell_com_df_gp.selected is True
    ].copy()

    # Metabolite
    adata.uns[ccc_key]["cell_com_df_m"]["selected"] = (
        (FDR_values_m < threshold) & (C_values_m > 0)
        if test == "non-parametric"
        else (FDR_values_m < threshold)
    )
    cell_com_df_m = adata.uns[ccc_key]["cell_com_df_m"]
    adata.uns[ccc_key]["cell_com_df_m_sig"] = cell_com_df_m[cell_com_df_m.selected is True].copy()

    return


def compute_interacting_cell_scores(
    adata: str | AnnData,
    center_counts_for_np_test: bool | None = False,
    test: Literal["parametric"] | Literal["non-parametric"] | Literal["both"] | None = "both",
    restrict_significance: Literal["gene pairs"]
    | Literal["metabolites"]
    | Literal["both"]
    | None = "both",
    compute_significance: Literal["parametric"]
    | Literal["non-parametric"]
    | Literal["both"]
    | None = "both",
    M: int | None = 1000,
    seed: int | None = 42,
    check_analytic_null: bool | None = False,
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: bool | None = False,
):
    """
    Compute interacting cell scores for gene pairs and metabolites.

    Parameters
    ----------
    adata : AnnData or str
        AnnData object containing:
        - ``uns['model']`` and ``uns['mean']`` (expression normalization model)
        - ``uns['gene_pairs']``, ``uns['gene_pairs_per_metabolite']``
        - ``obsp['weights']``: sparse spatial weight matrix
        - ``uns['ccc_results']`` (for significance filtering)
    center_counts_for_np_test : bool, optional
        If True, center/normalize counts prior to the non-parametric test.
    test : {"parametric", "non-parametric", "both"}
        Which interacting cell score tests to compute.
    restrict_significance : {"gene pairs", "metabolites", "both"}
        Use only significant gene pairs/metabolites from CCC results.
    compute_significance : {"parametric", "non-parametric", "both"}
        Whether to compute significance (p-values, FDR) in each test.
    M : int, default 1000
        Number of permutations for the non-parametric test.
    seed : int, default 42
        Random seed for permutation reproducibility.
    check_analytic_null : bool, default False
        If True, evaluate the analytic null distribution during permutations.
    device : torch.device
        CPU or GPU device for tensor operations.
    verbose : bool, default False
        Print status updates.

    Returns
    -------
    None
        Results are stored in ``adata.uns['interacting_cell_results']``.
    """
    start = time.time()
    if verbose:
        print("Computing gene pair and metabolite scores...")

    adata.uns["interacting_cell_results"] = {}

    model = adata.uns["model"]
    mean = adata.uns["mean"]

    if test not in ["both", "parametric", "non-parametric"]:
        raise ValueError(
            'The "test" variable should be one of ["both", "parametric", "non-parametric"].'
        )

    if restrict_significance is not None and restrict_significance not in [
        "both",
        "gene pairs",
        "metabolites",
    ]:
        raise ValueError(
            'The "restrict_significance" variable should be one of ["both", "gene pairs", "metabolites"].'
        )

    if compute_significance is not None and compute_significance not in [
        "both",
        "parametric",
        "non-parametric",
    ]:
        raise ValueError(
            'The "compute_significance" variable should be one of ["both", "parametric", "non-parametric"].'
        )

    sample_specific = "sample_key" in adata.uns

    layer_key_p_test = adata.uns.get("layer_key_p_test", None)
    layer_key_np_test = adata.uns.get("layer_key_np_test", None)
    use_raw = (layer_key_p_test == "use_raw") and (layer_key_np_test == "use_raw")

    gene_pairs = adata.uns.get("gene_pairs", None)
    gene_pairs_per_metabolite = adata.uns["gene_pairs_per_metabolite"]

    def to_tuple(x):
        # Recursively convert lists to tuples
        if isinstance(x, list):
            return tuple(to_tuple(i) for i in x)
        return x

    metabolite_gene_pair_df = pd.DataFrame.from_dict(
        gene_pairs_per_metabolite, orient="index"
    ).reset_index()
    metabolite_gene_pair_df = metabolite_gene_pair_df.rename(columns={"index": "metabolite"})
    metabolite_gene_pair_df["gene_pair"] = metabolite_gene_pair_df["gene_pair"].apply(
        lambda arr: [(to_tuple(gp[0]), to_tuple(gp[1])) for gp in arr]
    )
    metabolite_gene_pair_df["gene_type"] = metabolite_gene_pair_df["gene_type"].apply(
        lambda arr: [(to_tuple(gt[0]), to_tuple(gt[1])) for gt in arr]
    )
    metabolite_gene_pair_df = pd.concat(
        [
            metabolite_gene_pair_df["metabolite"],
            metabolite_gene_pair_df.explode("gene_pair")["gene_pair"],
            metabolite_gene_pair_df.explode("gene_type")["gene_type"],
        ],
        axis=1,
    ).reset_index(drop=True)

    if "LR_database" in adata.uns:
        LR_database = adata.uns["LR_database"]
        df_merged = pd.merge(
            metabolite_gene_pair_df,
            LR_database,
            left_on="metabolite",
            right_on="interaction_name",
            how="left",
        )
        LR_df = df_merged.dropna(subset=["pathway_name"])
        metabolite_gene_pair_df["metabolite"][
            metabolite_gene_pair_df.metabolite.isin(LR_df.metabolite)
        ] = LR_df["pathway_name"]

    if restrict_significance in ["both", "gene pairs"]:
        cell_com_gp_df = adata.uns["ccc_results"]["cell_com_df_gp_sig"].copy()
        cell_com_gp_df[["Gene 1", "Gene 2"]] = cell_com_gp_df[["Gene 1", "Gene 2"]].map(
            lambda x: tuple(x) if isinstance(x, list) else x
        )

        gene_pairs_set = {tuple(x) for x in cell_com_gp_df[["Gene 1", "Gene 2"]].values}
        metabolite_gene_pair_df = metabolite_gene_pair_df[
            metabolite_gene_pair_df["gene_pair"].isin(gene_pairs_set)
        ]

    if restrict_significance in ["both", "metabolites"]:
        cell_com_m_df = adata.uns["ccc_results"]["cell_com_df_m_sig"].copy()
        metabolite_set = set(cell_com_m_df["Metabolite"].values)
        metabolite_gene_pair_df = metabolite_gene_pair_df[
            metabolite_gene_pair_df["metabolite"].isin(metabolite_set)
        ]

    genes = adata.uns["genes"]
    gene_pairs_sig = []
    if gene_pairs:
        for g1, g2 in gene_pairs:
            g1 = tuple(g1) if isinstance(g1, list) else g1
            g2 = tuple(g2) if isinstance(g2, list) else g2
            if not metabolite_gene_pair_df[metabolite_gene_pair_df["gene_pair"] == (g1, g2)].empty:
                gene_pairs_sig.append((g1, g2))

    adata.uns["gene_pairs_sig"] = gene_pairs_sig

    gene_pairs_sig_ind = []
    for g1, g2 in gene_pairs_sig:
        idx1 = tuple([genes.index(g) for g in g1]) if isinstance(g1, tuple) else genes.index(g1)
        idx2 = tuple([genes.index(g) for g in g2]) if isinstance(g2, tuple) else genes.index(g2)
        gene_pairs_sig_ind.append((idx1, idx2))

    adata.uns["gene_pairs_sig_ind"] = gene_pairs_sig_ind

    if "barcode_key" in adata.uns:
        barcode_key = adata.uns["barcode_key"]
        cells = pd.Series(adata.obs[barcode_key].tolist())
    else:
        cells = adata.obs_names if not use_raw else adata.raw.obs_names

    # Compute weights
    weights = make_weights_non_redundant(adata.obsp["weights"]).tocoo()
    weights = torch.sparse_coo_tensor(
        torch.tensor(np.vstack((weights.row, weights.col)), dtype=torch.long, device=device),
        torch.tensor(weights.data, dtype=torch.float64, device=device),
        torch.Size(weights.shape),
        device=device,
    )

    gene_pair_dict = {}
    for metabolite, group in metabolite_gene_pair_df.groupby("metabolite"):
        idxs = (
            group["gene_pair"]
            .apply(lambda gp: gene_pairs_sig.index(gp) if gp in gene_pairs_sig else None)
            .dropna()
            .tolist()
        )
        idxs = [int(ind) for ind in idxs if ind is not None]
        if idxs:
            gene_pair_dict[metabolite] = idxs
    metabolites = list(gene_pair_dict.keys())

    adata.uns["metabolites"] = metabolites

    gene_pairs_sig_names = [
        "_".join("_".join(g) if isinstance(g, tuple) else g for g in gp) for gp in gene_pairs_sig
    ]

    adata.uns["gene_pairs_sig_names"] = gene_pairs_sig_names

    if test in ["parametric", "both"]:
        if verbose:
            print("Running the parametric test...")

        adata.uns["interacting_cell_results"]["p"] = {"gp": {}, "m": {}}

        Wtot2 = torch.tensor((weights.data**2).sum(), device=device)

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_sig_ind:
            if isinstance(idx1, tuple):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, tuple):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        _lazy_import_hotspot()
        counts_1 = standardize_counts(adata, counts_1, model, num_umi, sample_specific)
        _lazy_import_hotspot()
        counts_2 = standardize_counts(adata, counts_2, model, num_umi, sample_specific)

        # Compute CCC scores
        WX2t = torch.sparse.mm(weights, counts_2.T)
        WtX2t = torch.sparse.mm(weights.transpose(0, 1), counts_2.T)
        cs_gp = (counts_1.T * WX2t) + (counts_1.T * WtX2t)
        same_gene_mask = torch.tensor([g1 == g2 for g1, g2 in gene_pairs_sig], device=device)
        cs_gp[:, same_gene_mask] = cs_gp[:, same_gene_mask] / 2
        adata.uns["interacting_cell_results"]["p"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()

        # Compute metabolite-level scores
        cs_m = compute_metabolite_cs(cs_gp, gene_pair_dict, interacting_cell_scores=True)
        adata.uns["interacting_cell_results"]["p"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        if compute_significance in ["parametric", "both"]:
            # Compute second moments
            WX1t = torch.sparse.mm(weights, counts_1.T)
            WtX1t = torch.sparse.mm(weights.transpose(0, 1), counts_1.T)
            eg2_a = (WX1t + WtX1t).pow(2)
            eg2_b = (WX2t + WtX2t).pow(2)
            eg2s_gp = (eg2_a, eg2_b)

            Z_gp, Z_m = compute_p_int_cell_results_no_ct(
                cs_gp, cs_m, gene_pairs_sig_ind, Wtot2, eg2s_gp, gene_pair_dict
            )

            Z_gp_np = Z_gp.detach().cpu().numpy()
            Z_m_np = Z_m.detach().cpu().numpy()
            # Compute p-values and FDRs
            Z_pvals_gp = norm.sf(Z_gp_np)
            Z_pvals_m = norm.sf(Z_m_np)
            FDR_gp = multipletests(Z_pvals_gp.flatten(), method="fdr_bh")[1].reshape(
                Z_pvals_gp.shape
            )
            FDR_m = multipletests(Z_pvals_m.flatten(), method="fdr_bh")[1].reshape(Z_pvals_m.shape)

            adata.uns["interacting_cell_results"]["p"]["gp"]["Z"] = Z_gp_np
            adata.uns["interacting_cell_results"]["p"]["gp"]["Z_pval"] = Z_pvals_gp
            adata.uns["interacting_cell_results"]["p"]["gp"]["Z_FDR"] = FDR_gp
            adata.uns["interacting_cell_results"]["p"]["m"]["Z"] = Z_m_np
            adata.uns["interacting_cell_results"]["p"]["m"]["Z_pval"] = Z_pvals_m
            adata.uns["interacting_cell_results"]["p"]["m"]["Z_FDR"] = FDR_m

            # P-value
            mask_gp = adata.uns["interacting_cell_results"]["p"]["gp"]["Z_pval"] < 0.05
            mask_m = adata.uns["interacting_cell_results"]["p"]["m"]["Z_pval"] < 0.05

            cs_gp_sig = adata.uns["interacting_cell_results"]["p"]["gp"]["cs"].copy()
            cs_m_sig = adata.uns["interacting_cell_results"]["p"]["m"]["cs"].copy()

            cs_gp_sig[~mask_gp] = np.nan
            cs_m_sig[~mask_m] = np.nan
            adata.uns["interacting_cell_results"]["p"]["gp"]["cs_sig_pval"] = cs_gp_sig
            adata.uns["interacting_cell_results"]["p"]["m"]["cs_sig_pval"] = cs_m_sig

            # FDR
            mask_gp = adata.uns["interacting_cell_results"]["p"]["gp"]["Z_FDR"] < 0.05
            mask_m = adata.uns["interacting_cell_results"]["p"]["m"]["Z_FDR"] < 0.05

            cs_gp_sig = adata.uns["interacting_cell_results"]["p"]["gp"]["cs"].copy()
            cs_m_sig = adata.uns["interacting_cell_results"]["p"]["m"]["cs"].copy()

            cs_gp_sig[~mask_gp] = np.nan
            cs_m_sig[~mask_m] = np.nan
            adata.uns["interacting_cell_results"]["p"]["gp"]["cs_sig_FDR"] = cs_gp_sig
            adata.uns["interacting_cell_results"]["p"]["m"]["cs_sig_FDR"] = cs_m_sig

        if verbose:
            print("Parametric test finished.")

    if test in ["non-parametric", "both"]:
        if verbose:
            print("Running the non-parametric test...")

        adata.uns["interacting_cell_results"]["np"] = {"gp": {}, "m": {}}

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_sig_ind:
            if isinstance(idx1, tuple):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, tuple):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        if center_counts_for_np_test:
            num_umi = counts.sum(dim=0)
            _lazy_import_hotspot()
            counts_1 = standardize_counts(adata, counts_1, model, num_umi, sample_specific)
            _lazy_import_hotspot()
            counts_2 = standardize_counts(adata, counts_2, model, num_umi, sample_specific)

        n_cells = counts_1.shape[1]
        same_gene_mask = torch.tensor([g1 == g2 for g1, g2 in gene_pairs_sig], device=device)

        if center_counts_for_np_test and test == "both":
            adata.uns["interacting_cell_results"]["np"]["gp"]["cs"] = np.array(
                adata.uns["interacting_cell_results"]["p"]["gp"]["cs"]
            )
            adata.uns["interacting_cell_results"]["np"]["m"]["cs"] = np.array(
                adata.uns["interacting_cell_results"]["p"]["m"]["cs"]
            )
        else:
            WX2t = torch.sparse.mm(weights, counts_2.T)
            WtX2t = torch.sparse.mm(weights.transpose(0, 1), counts_2.T)
            cs_gp = (counts_1.T * WX2t) + (counts_1.T * WtX2t)
            cs_gp[:, same_gene_mask] = cs_gp[:, same_gene_mask] / 2
            adata.uns["interacting_cell_results"]["np"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()
            cs_m = compute_metabolite_cs(cs_gp, gene_pair_dict, interacting_cell_scores=True)
            adata.uns["interacting_cell_results"]["np"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        if compute_significance in ["non-parametric", "both"]:
            perm_cs_gp_a = torch.zeros(
                (n_cells, counts_1.shape[0], M), dtype=torch.float64, device=device
            )
            perm_cs_gp_b = torch.zeros_like(perm_cs_gp_a)
            perm_cs_m_a = torch.zeros(
                (n_cells, len(gene_pair_dict), M), dtype=torch.float64, device=device
            )
            perm_cs_m_b = torch.zeros_like(perm_cs_m_a)

            if check_analytic_null:
                gp_zs_perm_array = torch.zeros_like(perm_cs_gp_a)
                gp_pvals_perm_array = torch.zeros_like(perm_cs_gp_a)
                m_zs_perm_array = torch.zeros_like(perm_cs_m_a)
                m_pvals_perm_array = torch.zeros_like(perm_cs_m_a)

            torch.manual_seed(seed)
            for i in tqdm(range(M), desc="Permutation test"):
                idx = torch.randperm(n_cells, device=device)

                c1_perm_a = counts_1.clone()
                c2_perm_a = counts_2[:, idx]
                c1_perm_a[same_gene_mask] = counts_1[same_gene_mask, :][:, idx]

                WX2t_a = torch.sparse.mm(weights, c2_perm_a.T)
                WtX2t_a = torch.sparse.mm(weights.transpose(0, 1), c2_perm_a.T)
                cs_a = (c1_perm_a.T * WX2t_a) + (c1_perm_a.T * WtX2t_a)
                cs_a[:, same_gene_mask] = cs_a[:, same_gene_mask] / 2
                perm_cs_gp_a[:, :, i] = cs_a

                cs_m_a = compute_metabolite_cs(cs_a, gene_pair_dict, interacting_cell_scores=True)
                perm_cs_m_a[:, :, i] = cs_m_a

                c2_perm_b = counts_2.clone()
                c1_perm_b = counts_1[:, idx]
                c2_perm_b[same_gene_mask] = counts_2[same_gene_mask, :][:, idx]

                WX2t_b = torch.sparse.mm(weights, c2_perm_b.T)
                WtX2t_b = torch.sparse.mm(weights.transpose(0, 1), c2_perm_b.T)
                cs_b = (c1_perm_b.T * WX2t_b) + (c1_perm_b.T * WtX2t_b)
                cs_b[:, same_gene_mask] = cs_b[:, same_gene_mask] / 2
                perm_cs_gp_b[:, :, i] = cs_b

                cs_m_b = compute_metabolite_cs(cs_b, gene_pair_dict, interacting_cell_scores=True)
                perm_cs_m_b[:, :, i] = cs_m_b

                if check_analytic_null:
                    Z_gp_perm, Z_m_perm = compute_p_results(
                        (cs_a, cs_b),
                        (cs_m_a, cs_m_b),
                        gene_pairs_ind,
                        Wtot2,
                        eg2s_gp,
                        gene_pair_dict,
                    )
                    gp_zs_perm_array[:, :, i] = Z_gp_perm
                    gp_pvals_perm_array[:, :, i] = torch.tensor(
                        norm.sf(Z_gp_perm.cpu().numpy()), device=device
                    )
                    m_zs_perm_array[:, :, i] = Z_m_perm
                    m_pvals_perm_array[:, :, i] = torch.tensor(
                        norm.sf(Z_m_perm.cpu().numpy()), device=device
                    )

            adata.uns["interacting_cell_results"]["np"]["gp"]["perm_cs_a"] = (
                perm_cs_gp_a.detach().cpu().numpy()
            )
            adata.uns["interacting_cell_results"]["np"]["gp"]["perm_cs_b"] = (
                perm_cs_gp_b.detach().cpu().numpy()
            )
            adata.uns["interacting_cell_results"]["np"]["m"]["perm_cs_a"] = (
                perm_cs_m_a.detach().cpu().numpy()
            )
            adata.uns["interacting_cell_results"]["np"]["m"]["perm_cs_b"] = (
                perm_cs_m_b.detach().cpu().numpy()
            )

            x_gp_a = (perm_cs_gp_a > cs_gp[:, :, None]).sum(dim=2)
            x_gp_b = (perm_cs_gp_b > cs_gp[:, :, None]).sum(dim=2)
            x_m_a = (perm_cs_m_a > cs_m[:, :, None]).sum(dim=2)
            x_m_b = (perm_cs_m_b > cs_m[:, :, None]).sum(dim=2)

            pvals_gp_a = (x_gp_a + 1).float() / (M + 1)
            pvals_gp_b = (x_gp_b + 1).float() / (M + 1)
            pvals_m_a = (x_m_a + 1).float() / (M + 1)
            pvals_m_b = (x_m_b + 1).float() / (M + 1)

            pvals_gp = torch.where(pvals_gp_a > pvals_gp_b, pvals_gp_a, pvals_gp_b)
            pvals_m = torch.where(pvals_m_a > pvals_m_b, pvals_m_a, pvals_m_b)

            pvals_gp = pvals_gp.cpu().numpy()
            pvals_m = pvals_m.cpu().numpy()

            adata.uns["interacting_cell_results"]["np"]["gp"]["pval"] = pvals_gp
            adata.uns["interacting_cell_results"]["np"]["gp"]["FDR"] = multipletests(
                pvals_gp.flatten(), method="fdr_bh"
            )[1].reshape(pvals_gp.shape)
            adata.uns["interacting_cell_results"]["np"]["m"]["pval"] = pvals_m
            adata.uns["interacting_cell_results"]["np"]["m"]["FDR"] = multipletests(
                pvals_m.flatten(), method="fdr_bh"
            )[1].reshape(pvals_m.shape)

            if check_analytic_null:
                adata.uns["interacting_cell_results"]["np"]["analytic_null"] = {
                    "gp_zs_perm": gp_zs_perm_array.detach().cpu().numpy(),
                    "gp_pvals_perm": gp_pvals_perm_array.detach().cpu().numpy(),
                    "m_zs_perm": m_zs_perm_array.detach().cpu().numpy(),
                    "m_pvals_perm": m_pvals_perm_array.detach().cpu().numpy(),
                }

            # P-value
            mask_gp = adata.uns["interacting_cell_results"]["np"]["gp"]["pval"] < 0.05
            mask_m = adata.uns["interacting_cell_results"]["np"]["m"]["pval"] < 0.05

            cs_gp_sig = adata.uns["interacting_cell_results"]["np"]["gp"]["cs"].copy()
            cs_m_sig = adata.uns["interacting_cell_results"]["np"]["m"]["cs"].copy()

            cs_gp_sig[~mask_gp] = np.nan
            cs_m_sig[~mask_m] = np.nan
            adata.uns["interacting_cell_results"]["np"]["gp"]["cs_sig_pval"] = cs_gp_sig
            adata.uns["interacting_cell_results"]["np"]["m"]["cs_sig_pval"] = cs_m_sig

            # FDR
            mask_gp = adata.uns["interacting_cell_results"]["np"]["gp"]["FDR"] < 0.05
            mask_m = adata.uns["interacting_cell_results"]["np"]["m"]["FDR"] < 0.05

            cs_gp_sig = adata.uns["interacting_cell_results"]["np"]["gp"]["cs"].copy()
            cs_m_sig = adata.uns["interacting_cell_results"]["np"]["m"]["cs"].copy()

            cs_gp_sig[~mask_gp] = np.nan
            cs_m_sig[~mask_m] = np.nan
            adata.uns["interacting_cell_results"]["np"]["gp"]["cs_sig_FDR"] = cs_gp_sig
            adata.uns["interacting_cell_results"]["np"]["m"]["cs_sig_FDR"] = cs_m_sig

        if verbose:
            print("Non-parametric test finished.")

    if verbose:
        print(
            "Finished computing gene pair and metabolite scores in %.3f seconds"
            % (time.time() - start)
        )

    return


def compute_ct_interacting_cell_scores(
    adata: str | AnnData,
    center_counts_for_np_test: bool | None = False,
    test: Literal["parametric"] | Literal["non-parametric"] | Literal["both"] | None = "both",
    restrict_significance: Literal["gene pairs"]
    | Literal["metabolites"]
    | Literal["both"]
    | None = "both",
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    verbose: bool | None = False,
):
    """
    Compute cell-type–aware interacting cell scores for gene pairs and metabolites.

    Parameters
    ----------
    adata : AnnData or str
        Must contain:
        - ``uns['model']``, ``uns['mean']``
        - ``uns['cell_type_key']`` and ``obs[cell_type_key]`` for cell types
        - ``uns['gene_pairs']``, ``uns['gene_pairs_per_ct_pair']``
        - ``uns['gene_pairs_per_metabolite']``
        - ``uns['ct_ccc_results']`` with significance information
        - ``obsp['weights']`` (spatial proximity matrix)
    center_counts_for_np_test : bool, default False
        Whether to standardize counts before the non-parametric test.
    test : {"parametric", "non-parametric", "both"}
        Which statistical test(s) to run.
    restrict_significance : {"gene pairs", "metabolites", "both"}
        Only use cell-type-pair interactions that were significant in cell-type-aware CCC results.
    device : torch.device
        CPU or GPU device for PyTorch computations.
    verbose : bool, default False
        Print detailed progress messages.
    """
    start = time.time()
    if verbose:
        print("Computing cell type-aware gene pair and metabolite scores...")

    adata.uns["ct_interacting_cell_results"] = {}

    model = adata.uns["model"]
    mean = adata.uns["mean"]

    if test not in ["both", "parametric", "non-parametric"]:
        raise ValueError(
            'The "test" variable should be one of ["both", "parametric", "non-parametric"].'
        )

    if restrict_significance not in ["both", "gene pairs", "metabolites"]:
        raise ValueError(
            'The "restrict_significance" variable should be one of ["both", "gene pairs", "metabolites"].'
        )

    sample_specific = "sample_key" in adata.uns

    layer_key_p_test = adata.uns.get("layer_key_p_test", None)
    layer_key_np_test = adata.uns.get("layer_key_np_test", None)
    use_raw = (layer_key_p_test == "use_raw") and (layer_key_np_test == "use_raw")

    obs = adata.raw.obs if use_raw else adata.obs
    cells = (
        adata.raw.obs.index.values.astype(str) if use_raw else adata.obs_names.values.astype(str)
    )

    gene_pairs = adata.uns.get("gene_pairs", None)
    gene_pairs_per_ct_pair = adata.uns.get("gene_pairs_per_ct_pair", None)

    gp_metab = adata.uns["gene_pairs_per_metabolite"]
    metabolite_gene_pair_df = (
        pd.DataFrame.from_dict(gp_metab, orient="index")
        .rename_axis("metabolite")
        .explode(["gene_pair", "gene_type"])
        .reset_index()
    )

    if "LR_database" in adata.uns:
        merged = metabolite_gene_pair_df.merge(
            adata.uns["LR_database"], left_on="metabolite", right_on="interaction_name", how="left"
        )
        LR_df = merged.dropna(subset=["pathway_name"])
        metabolite_gene_pair_df.loc[
            metabolite_gene_pair_df.metabolite.isin(LR_df.metabolite), "metabolite"
        ] = LR_df["pathway_name"].values

    cell_type_pairs = adata.uns.get("cell_type_pairs")
    cell_type_pairs = [tuple(x) for x in cell_type_pairs]

    cell_com_gp_df = adata.uns["ct_ccc_results"]["cell_com_df_gp_sig"].copy()
    if restrict_significance in ["both", "gene pairs"]:
        ct_pairs_gp_set = {tuple(x) for x in cell_com_gp_df[["Cell Type 1", "Cell Type 2"]].values}
        cell_type_pairs = [ct_pair for ct_pair in cell_type_pairs if ct_pair in ct_pairs_gp_set]

        cell_com_gp_df[["Gene 1", "Gene 2"]] = cell_com_gp_df[["Gene 1", "Gene 2"]].map(
            lambda x: tuple(x) if isinstance(x, list) else x
        )

        gene_pairs_set = {tuple(x) for x in cell_com_gp_df[["Gene 1", "Gene 2"]].values}
        metabolite_gene_pair_df = metabolite_gene_pair_df[
            metabolite_gene_pair_df["gene_pair"].isin(gene_pairs_set)
        ]

    cell_com_m_df = adata.uns["ct_ccc_results"]["cell_com_df_m_sig"].copy()
    if restrict_significance in ["both", "metabolites"]:
        ct_pairs_m_set = {tuple(x) for x in cell_com_m_df[["Cell Type 1", "Cell Type 2"]].values}
        missing_ct_pairs = [
            ct_pair for ct_pair in ct_pairs_m_set if ct_pair not in cell_type_pairs
        ]
        if len(missing_ct_pairs) > 0:
            warnings.warn(
                f'The following cell type pairs are not included in the "cell_type_pairs" set: {missing_ct_pairs}',
                stacklevel=2,
            )

        metabolite_set = set(cell_com_m_df["metabolite"].values)
        metabolite_gene_pair_df = metabolite_gene_pair_df[
            metabolite_gene_pair_df["metabolite"].isin(metabolite_set)
        ]

    if metabolite_gene_pair_df.empty:
        if restrict_significance == "both":
            raise ValueError(
                "There are no significant gene pairs that belong to a significant metabolite."
            )
        if restrict_significance == "gene pairs":
            raise ValueError("There are no significant gene pairs.")
        if restrict_significance == "metabolites":
            raise ValueError("There are no significant metabolites.")

    genes = adata.uns["genes"]
    gene_pairs_sig = []
    if gene_pairs:
        for g1, g2 in gene_pairs:
            g1 = tuple(g1) if isinstance(g1, list) else g1
            g2 = tuple(g2) if isinstance(g2, list) else g2
            if not metabolite_gene_pair_df[metabolite_gene_pair_df["gene_pair"] == (g1, g2)].empty:
                gene_pairs_sig.append((g1, g2))

    gene_pairs_sig_ind = []
    for pair in gene_pairs_sig:
        idx1 = (
            [genes.index(g) for g in pair[0]]
            if isinstance(pair[0], list)
            else genes.index(pair[0])
        )
        idx2 = (
            [genes.index(g) for g in pair[1]]
            if isinstance(pair[1], list)
            else genes.index(pair[1])
        )
        gene_pairs_sig_ind.append((idx1, idx2))

    cell_type_key = adata.uns.get("cell_type_key")
    cell_types = obs[cell_type_key]
    gene_pairs_per_ct_pair = adata.uns.get("gene_pairs_per_ct_pair", {})

    weights = adata.obsp["weights"]

    used_ct_pairs = list({ct for cell_type_pair in cell_type_pairs for ct in cell_type_pair})
    all_cell_types = set(cell_types.unique())
    used_ct_pairs_set = set(used_ct_pairs)
    if used_ct_pairs_set < all_cell_types:
        keep_mask = cell_types[cells].isin(used_ct_pairs).values
        keep_indices = np.where(keep_mask)[0]
        weights = weights[keep_indices][:, keep_indices]
        cells = cells[keep_indices]
        cell_types = cell_types.loc[
            cells
        ]  # Eventually only keep the cell type pairs with at least one significant gene pair

    weights_ct_pairs = create_weights_ct_pairs(
        weights.tocoo(), cell_types, cell_type_pairs, device
    )

    gene_pairs_per_ct_pair_sig = {}
    for ct_pair in gene_pairs_per_ct_pair.keys():
        if ct_pair not in cell_type_pairs:
            continue
        cell_com_df_ct_pair = cell_com_gp_df[
            (cell_com_gp_df["Cell Type 1"] == ct_pair[0])
            & (cell_com_gp_df["Cell Type 2"] == ct_pair[1])
        ]
        gene_pairs_per_ct_pair_sig[ct_pair] = [
            tuple(x) for x in cell_com_df_ct_pair[["Gene 1", "Gene 2"]].values
        ]

    # Cell-type pair-specific indices
    gene_pairs_ind_per_ct_pair_sig = defaultdict(list)
    gene_pairs_per_ct_pair_sig_ind = defaultdict(list)
    for ct_pair, gpairs in gene_pairs_per_ct_pair_sig.items():
        for pair in gpairs:
            if pair not in gene_pairs_sig:
                continue
            idx = gene_pairs_sig.index(pair)
            gene_pairs_ind_per_ct_pair_sig[ct_pair].append(gene_pairs_sig_ind[idx])
            gene_pairs_per_ct_pair_sig_ind[ct_pair].append(idx)

    adata.uns["gene_pairs_ind_per_ct_pair_sig"] = dict(gene_pairs_ind_per_ct_pair_sig)
    adata.uns["gene_pairs_per_ct_pair_sig_ind"] = dict(gene_pairs_per_ct_pair_sig_ind)

    def make_hashable(pair):
        return tuple(tuple(x) if isinstance(x, list) else x for x in pair)

    gene_pairs_sig_ind_set = {make_hashable(pair) for pair in gene_pairs_sig_ind}
    ct_specific_gene_pairs = [
        i
        for i, pairs in enumerate(gene_pairs_ind_per_ct_pair_sig.values())
        if {make_hashable(pair) for pair in pairs} < gene_pairs_sig_ind_set
    ]

    gene_pair_dict = {}
    for metabolite, group in metabolite_gene_pair_df.groupby("metabolite"):
        idxs = (
            group["gene_pair"]
            .apply(lambda gp: gene_pairs_sig.index(gp) if gp in gene_pairs_sig else None)
            .dropna()
            .tolist()
        )
        idxs = [int(ind) for ind in idxs if ind is not None]
        if idxs:
            gene_pair_dict[metabolite] = idxs
    metabolites = list(gene_pair_dict.keys())

    gene_pair_to_metabolite_indices = defaultdict(set)
    for met_idx, met in enumerate(metabolites):
        for gene_pair_idx in gene_pair_dict.get(met, []):
            gene_pair_to_metabolite_indices[gene_pair_idx].add(met_idx)

    ct_pair_to_metabolite_indices = {}
    for ct_pair, gene_pair_indices in gene_pairs_per_ct_pair_sig_ind.items():
        met_indices = set()
        for gi in gene_pair_indices:
            met_indices.update(gene_pair_to_metabolite_indices.get(gi, []))
        ct_pair_to_metabolite_indices[ct_pair] = sorted(met_indices)

    def concat_tuple_elements(t, sep="_"):
        return sep.join(
            s for item in t for s in (item if isinstance(item, list | tuple) else [item])
        )

    if test in ["parametric", "both"]:
        if verbose:
            print("Running the parametric test...")

        adata.uns["ct_interacting_cell_results"]["p"] = {"gp": {}, "m": {}}

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_sig_ind:
            if isinstance(idx1, list):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, list):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        counts_1 = standardize_ct_counts(
            adata[cells, :], counts_1, model, num_umi, sample_specific, cell_types
        )
        counts_2 = standardize_ct_counts(
            adata[cells, :], counts_2, model, num_umi, sample_specific, cell_types
        )

        cs_gp = torch.zeros(
            (len(cell_type_pairs), counts_1.shape[1], counts_1.shape[0]), device=counts_1.device
        )
        for ct_pair in range(len(cell_type_pairs)):
            W = weights_ct_pairs[ct_pair].coalesce()
            WX2t = torch.sparse.mm(W, counts_2.T)
            cs_gp[ct_pair] = counts_1.T * WX2t
        adata.uns["ct_interacting_cell_results"]["p"]["gp"]["cs"] = cs_gp.detach().cpu().numpy()

        cs_m = compute_metabolite_cs_ct(
            cs_gp,
            cell_type_key,
            gene_pair_dict,
            gene_pairs_per_ct_pair_sig_ind,
            ct_specific_gene_pairs,
            interacting_cell_scores=True,
        )
        adata.uns["ct_interacting_cell_results"]["p"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        column_names = []
        scores = []
        for i, ct_pair in enumerate(gene_pairs_per_ct_pair_sig_ind.keys()):
            gp_list = gene_pairs_per_ct_pair_sig_ind[ct_pair]
            for gp in gp_list:
                if np.all(
                    adata.uns["ct_interacting_cell_results"]["p"]["gp"]["cs"][i, :, gp] == 0
                ):
                    continue
                column_names.append(
                    f"{' - '.join(ct_pair)}: {concat_tuple_elements(gene_pairs_sig[gp])}"
                )
                scores.append(adata.uns["ct_interacting_cell_results"]["p"]["gp"]["cs"][i, :, gp])

        cs_gp_df = pd.DataFrame(
            {column_names[i]: array for i, array in enumerate(scores)}, index=cells
        )
        if used_ct_pairs_set < all_cell_types:
            cs_gp_df = cs_gp_df.reindex(adata.obs_names, fill_value=0)
        adata.obsm["ct_interacting_cell_results_p_gp_cs_df"] = cs_gp_df

        column_names = []
        scores = []
        for i, ct_pair in enumerate(ct_pair_to_metabolite_indices.keys()):
            metab_list = ct_pair_to_metabolite_indices[ct_pair]
            for metab in metab_list:
                if np.all(
                    adata.uns["ct_interacting_cell_results"]["p"]["m"]["cs"][i, :, metab] == 0
                ):
                    continue
                column_names.append(f"{' - '.join(ct_pair)}: {metabolites[metab]}")
                scores.append(
                    adata.uns["ct_interacting_cell_results"]["p"]["m"]["cs"][i, :, metab]
                )

        cs_m_df = pd.DataFrame(
            {column_names[i]: array for i, array in enumerate(scores)}, index=cells
        )
        if used_ct_pairs_set < all_cell_types:
            cs_m_df = cs_m_df.reindex(adata.obs_names, fill_value=0)
        adata.obsm["ct_interacting_cell_results_p_m_cs_df"] = cs_m_df

        if verbose:
            print("Parametric test finished.")

    if test in ["non-parametric", "both"]:
        if verbose:
            print("Running the non-parametric test...")

        adata.uns["ct_interacting_cell_results"]["np"] = {"gp": {}, "m": {}}

        # Load counts
        counts = counts_from_anndata(adata[cells, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)

        # Prepare counts_1 and counts_2
        counts_1 = []
        counts_2 = []
        for idx1, idx2 in gene_pairs_sig_ind:
            if isinstance(idx1, tuple):
                c1 = (
                    counts[idx1, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx1, :] + 1e-8).mean(dim=0))
                )
            else:
                c1 = counts[idx1, :]
            if isinstance(idx2, tuple):
                c2 = (
                    counts[idx2, :].mean(dim=0)
                    if mean == "algebraic"
                    else torch.exp(torch.log(counts[idx2, :] + 1e-8).mean(dim=0))
                )
            else:
                c2 = counts[idx2, :]
            counts_1.append(c1)
            counts_2.append(c2)

        counts_1 = torch.stack(counts_1)
        counts_2 = torch.stack(counts_2)

        if center_counts_for_np_test:
            num_umi = counts.sum(dim=0)
            _lazy_import_hotspot()
            counts_1 = standardize_counts(adata, counts_1, model, num_umi, sample_specific)
            _lazy_import_hotspot()
            counts_2 = standardize_counts(adata, counts_2, model, num_umi, sample_specific)

        if center_counts_for_np_test and test == "both":
            adata.uns["ct_interacting_cell_results"]["np"]["gp"]["cs"] = np.array(
                adata.uns["ct_interacting_cell_results"]["p"]["gp"]["cs"]
            )
            adata.uns["ct_interacting_cell_results"]["np"]["m"]["cs"] = np.array(
                adata.uns["ct_interacting_cell_results"]["p"]["m"]["cs"]
            )
        else:
            cs_gp = torch.zeros(
                (len(cell_type_pairs), counts_1.shape[1], counts_1.shape[0]),
                device=counts_1.device,
            )
            for ct_pair in range(len(cell_type_pairs)):
                W = weights_ct_pairs[ct_pair].coalesce()
                WX2t = torch.sparse.mm(W, counts_2.T)
                cs_gp[ct_pair] = counts_1.T * WX2t
            adata.uns["ct_interacting_cell_results"]["np"]["gp"]["cs"] = (
                cs_gp.detach().cpu().numpy()
            )
            cs_m = compute_metabolite_cs_ct(
                cs_gp,
                cell_type_key,
                gene_pair_dict,
                gene_pairs_per_ct_pair_sig_ind,
                ct_specific_gene_pairs,
                interacting_cell_scores=True,
            )
            adata.uns["ct_interacting_cell_results"]["np"]["m"]["cs"] = cs_m.detach().cpu().numpy()

        column_names = []
        scores = []
        for i, ct_pair in enumerate(gene_pairs_per_ct_pair_sig_ind.keys()):
            gp_list = gene_pairs_per_ct_pair_sig_ind[ct_pair]
            for gp in gp_list:
                if np.all(
                    adata.uns["ct_interacting_cell_results"]["np"]["gp"]["cs"][i, :, gp] == 0
                ):
                    continue
                column_names.append(
                    f"{' - '.join(ct_pair)}: {concat_tuple_elements(gene_pairs_sig[gp])}"
                )
                scores.append(adata.uns["ct_interacting_cell_results"]["np"]["gp"]["cs"][i, :, gp])

        cs_gp_df = pd.DataFrame(
            {column_names[i]: array for i, array in enumerate(scores)}, index=cells
        )
        if used_ct_pairs_set < all_cell_types:
            cs_gp_df = cs_gp_df.reindex(adata.obs_names, fill_value=0)
        adata.obsm["ct_interacting_cell_results_np_gp_cs_df"] = cs_gp_df

        column_names = []
        scores = []
        for i, ct_pair in enumerate(ct_pair_to_metabolite_indices.keys()):
            metab_list = ct_pair_to_metabolite_indices[ct_pair]
            for metab in metab_list:
                if np.all(
                    adata.uns["ct_interacting_cell_results"]["np"]["m"]["cs"][i, :, metab] == 0
                ):
                    continue
                column_names.append(f"{' - '.join(ct_pair)}: {metabolites[metab]}")
                scores.append(
                    adata.uns["ct_interacting_cell_results"]["np"]["m"]["cs"][i, :, metab]
                )

        cs_m_df = pd.DataFrame(
            {column_names[i]: array for i, array in enumerate(scores)}, index=cells
        )
        if used_ct_pairs_set < all_cell_types:
            cs_m_df = cs_m_df.reindex(adata.obs_names, fill_value=0)
        adata.obsm["ct_interacting_cell_results_np_m_cs_df"] = cs_m_df

        if verbose:
            print("Non-parametric test finished.")

    if verbose:
        print(
            "Finished computing cell type-aware gene pair and metabolite scores in %.3f seconds"
            % (time.time() - start)
        )

    return


def compute_metabolite_cs(
    cs_gp: torch.Tensor, gene_pair_dict: dict, interacting_cell_scores: bool = False
) -> torch.Tensor:
    """
    Computes metabolite-level communication scores from gene-pair scores.

    Parameters
    ----------
    cs_gp : torch.Tensor
        - If interacting_cell_scores is False: shape (gene_pairs,)
        - If interacting_cell_scores is True: shape (cells, gene_pairs)
    gene_pair_dict : dict
        Maps metabolite names to a list of indices (ints) referring to gene-pairs.
    interacting_cell_scores : bool, optional
        Whether cs_gp contains per-cell scores.

    Returns
    -------
    cs_m : torch.Tensor
        - If interacting_cell_scores is False: shape (num_metabolites,)
        - If interacting_cell_scores is True: shape (cells, num_metabolites)
    """
    device = cs_gp.device
    scores = []

    for indices in gene_pair_dict.values():
        idx_tensor = torch.tensor(indices, device=device, dtype=torch.long)
        if interacting_cell_scores:
            summed = cs_gp[:, idx_tensor].sum(dim=1)  # shape: (cells,)
        else:
            summed = cs_gp[idx_tensor].sum()  # scalar
        scores.append(summed)

    if interacting_cell_scores:
        cs_m = torch.stack(scores, dim=1)  # shape: (cells, metabolites)
    else:
        cs_m = torch.stack(scores)  # shape: (metabolites,)

    return cs_m


def compute_metabolite_cs_ct(
    cs_gp,
    cell_type_key,
    gene_pair_dict,
    gene_pairs_per_ct_pair_ind=None,
    ct_specific_gene_pairs=None,
    interacting_cell_scores=False,
):
    if cell_type_key and ct_specific_gene_pairs:
        for i, ct_pair in enumerate(gene_pairs_per_ct_pair_ind.keys()):
            if i not in ct_specific_gene_pairs:
                continue
            mask_dim = 2 if interacting_cell_scores else 1
            mask = np.ones(cs_gp.shape[mask_dim], dtype=bool)
            mask[gene_pairs_per_ct_pair_ind[ct_pair]] = False
            if interacting_cell_scores:
                cs_gp[i, :, mask] = 0
            else:
                cs_gp[i, mask] = 0

    device = cs_gp.device
    scores = []

    for indices in gene_pair_dict.values():
        idx_tensor = torch.tensor(indices, device=device, dtype=torch.long)
        if interacting_cell_scores:
            summed = cs_gp[:, :, idx_tensor].sum(dim=2)  # shape: (cells,)
        else:
            summed = cs_gp[:, idx_tensor].sum(dim=1)  # scalar
        scores.append(summed)

    if interacting_cell_scores:
        cs_m = torch.stack(scores, dim=2)  # shape: (cells, metabolites)
    else:
        cs_m = torch.stack(scores, dim=1)

    return cs_m


def compute_metabolite_cs_old(
    cs_gp,
    cell_type_key,
    gene_pair_dict,
    gene_pairs_per_ct_pair_ind=None,
    ct_specific_gene_pairs=None,
    interacting_cell_scores=False,
):
    if cell_type_key and ct_specific_gene_pairs:
        for i, ct_pair in enumerate(gene_pairs_per_ct_pair_ind.keys()):
            if i not in ct_specific_gene_pairs:
                continue
            mask = np.ones(cs_gp.shape[1], dtype=bool)
            mask[gene_pairs_per_ct_pair_ind[ct_pair]] = False
            cs_gp[i, mask] = 0

    cells_metabolites = []
    for _metabolite, gene_pair_indices in gene_pair_dict.items():
        if interacting_cell_scores:
            summed_values = (
                cs_gp[:, :, gene_pair_indices].sum(axis=2)
                if cell_type_key
                else cs_gp[:, gene_pair_indices].sum(axis=1)
            )
            cells_metabolites.append(summed_values)
        else:
            summed_values = (
                cs_gp[:, gene_pair_indices].sum(axis=1)
                if cell_type_key
                else cs_gp[gene_pair_indices].sum(axis=0)
            )
            cells_metabolites.append(summed_values)
    if interacting_cell_scores:
        axis = 2 if cell_type_key else 1
    else:
        axis = 1 if cell_type_key else 0
    cs_m = np.stack(cells_metabolites, axis=axis)

    return cs_m


def ensure_tuple(x):
    return tuple(tuple(i) if isinstance(i, list) else i for i in x)


def compute_CCC_scores(
    counts_1: np.array,
    counts_2: np.array,
    weights: sparse.COO,
    gene_pairs: list,
):

    if len(weights.shape) == 3:
        scores = (counts_1.T * np.tensordot(weights, counts_2.T, axes=([2], [0]))).sum(axis=1)
    else:
        same_gene_mask = np.array([pair1 == pair2 for pair1, pair2 in gene_pairs])
        scores = (counts_1.T * (weights @ counts_2.T)).sum(axis=0) + (
            counts_1.T * (weights.T @ counts_2.T)
        ).sum(axis=0)
        scores[same_gene_mask] = scores[same_gene_mask] / 2

    return scores


def compute_int_CCC_scores(
    counts_1: np.array,
    counts_2: np.array,
    weights: sparse.COO,
    gene_pairs: list,
):

    if len(weights.shape) == 3:
        scores = counts_1.T * np.tensordot(weights, counts_2.T, axes=([2], [0]))
    else:
        same_gene_mask = np.array([pair1 == pair2 for pair1, pair2 in gene_pairs])
        scores = (counts_1.T * (weights @ counts_2.T)) + (counts_1.T * (weights.T @ counts_2.T))
        scores[:, same_gene_mask] = scores[:, same_gene_mask] / 2

    return scores


def get_ct_pair_weights(weights, cell_type_pairs, cell_types):

    w_nrow, w_ncol = weights.shape
    n_ct_pairs = len(cell_type_pairs)

    extract_weights_results = partial(
        extract_ct_pair_weights,
        weights=weights,
        cell_type_pairs=cell_type_pairs,
        cell_types=cell_types,
    )
    results = list(map(extract_weights_results, cell_type_pairs))

    w_new_data_all = [x[0] for x in results]
    w_new_coords_3d_all = [x[1] for x in results]
    w_new_coords_3d_all = np.hstack(w_new_coords_3d_all)
    w_new_data_all = np.concatenate(w_new_data_all)

    weights_ct_pairs = sparse.COO(
        w_new_coords_3d_all, w_new_data_all, shape=(n_ct_pairs, w_nrow, w_ncol)
    )

    return weights_ct_pairs


def extract_ct_pair_weights(ct_pair, weights, cell_type_pairs, cell_types):

    i = cell_type_pairs.index(ct_pair)

    ct_t, ct_u = cell_type_pairs[i]
    ct_t_mask = cell_types.values == ct_t
    ct_t_mask_coords = np.argwhere(ct_t_mask)
    ct_u_mask = cell_types.values == ct_u
    ct_u_mask_coords = np.argwhere(ct_u_mask)

    w_old_coords = weights.coords

    w_row_coords, w_col_coords = np.meshgrid(ct_t_mask_coords, ct_u_mask_coords, indexing="ij")
    w_row_coords = w_row_coords.ravel()
    w_col_coords = w_col_coords.ravel()
    w_new_coords = np.vstack((w_row_coords, w_col_coords))

    # w_matching_indices = np.where(np.all(np.isin(w_old_coords.T, w_new_coords.T), axis=1))[0]
    w_matching_indices = np.isin(w_old_coords[0], ct_t_mask_coords.flatten()) & np.isin(
        w_old_coords[1], ct_u_mask_coords.flatten()
    )
    w_new_data = weights.data[w_matching_indices]
    w_new_coords = w_old_coords[:, w_matching_indices]

    w_coord_3d = np.full(w_new_coords.shape[1], fill_value=i)
    w_new_coords_3d = np.vstack((w_coord_3d, w_new_coords))

    return (w_new_data, w_new_coords_3d)


def get_interacting_cell_type_pairs(x, weights, cell_types):
    ct_1, ct_2 = x

    ct_1_bin = cell_types == ct_1
    ct_2_bin = cell_types == ct_2

    weights = weights.tocsc()
    cell_types_weights = weights[ct_1_bin,][:, ct_2_bin]

    return bool(cell_types_weights.nnz)


def conditional_eg2_cellcom_gp(counts_1, counts_2, weights):

    if len(weights.shape) == 3:
        counts_1_sq = counts_1**2
        counts_2_sq = counts_2**2
        weights_sq_data = weights.data**2
        weights_sq = sparse.COO(weights.coords, weights_sq_data, shape=weights.shape)
        out_eg2_a = np.tensordot(counts_1_sq, weights_sq, axes=([1], [1])).sum(axis=2).T
        out_eg2_b = np.tensordot(weights_sq, counts_2_sq, axes=([2], [1])).sum(axis=1)
        out_eg2s = (out_eg2_a, out_eg2_b)
    else:
        out_eg2_a = (((weights + weights.T) @ counts_1.T) ** 2).sum(axis=0)
        out_eg2_b = (((weights + weights.T) @ counts_2.T) ** 2).sum(axis=0)
        out_eg2s = (out_eg2_a, out_eg2_b)

    return out_eg2s


def conditional_eg2_gp_score(counts, weights):

    counts_sq = counts**2
    if len(weights.shape) == 3:
        # weights_t = weights.transpose(axes=(0, 2, 1))
        # weights = weights + weights_t
        weights_sq_data = weights.data**2
        weights_sq = sparse.COO(weights.coords, weights_sq_data, shape=weights.shape)
        out_eg2_a = np.transpose(np.tensordot(counts_sq, weights_sq, axes=([1], [1])), (1, 2, 0))
        out_eg2_b = np.tensordot(weights_sq, counts_sq, axes=([2], [1]))
        out_eg2s = (out_eg2_a, out_eg2_b)
    else:
        out_eg2s = ((weights + weights.T) @ counts.T) ** 2

    return out_eg2s


def compute_ct_p_results(
    C_gp,
    C_m,
    gene_pairs_per_ct_pair_ind,
    ct_specific_gene_pairs,
    EG2_gp,
    cell_type_key,
    gene_pair_dict,
):

    EG2_gp = EG2_gp.unsqueeze(1).expand(-1, C_gp.shape[1]) if len(EG2_gp.shape) == 1 else EG2_gp

    stdG = torch.sqrt(EG2_gp)
    stdG[stdG == 0] = 1

    Z_gp = C_gp / stdG

    EG2_m = compute_metabolite_cs_ct(
        EG2_gp,
        cell_type_key,
        gene_pair_dict,
        gene_pairs_per_ct_pair_ind,
        ct_specific_gene_pairs,
        interacting_cell_scores=False,
    )
    if not isinstance(EG2_m, torch.Tensor):
        device = EG2_gp.device
        EG2_m = torch.tensor(EG2_m, device=device, dtype=torch.float64)

    stdG_m = torch.sqrt(EG2_m)
    stdG_m[stdG_m == 0] = 1

    Z_m = C_m / stdG_m

    return Z_gp, Z_m


def compute_p_results(C_gp, C_m, gene_pairs_ind, Wtot2, eg2s_gp, gene_pair_dict):

    device = Wtot2.device

    # Convert indices
    same_gene_mask = torch.tensor(
        [
            (isinstance(g1, int) and isinstance(g2, int) and g1 == g2)
            or (isinstance(g1, list) and isinstance(g2, list) and sorted(g1) == sorted(g2))
            for g1, g2 in gene_pairs_ind
        ],
        device=device,
    )

    # Unpack second moments
    EG2_a = eg2s_gp[0].clone()
    EG2_b = eg2s_gp[1].clone()
    EG2_a[same_gene_mask] = Wtot2
    EG2_b[same_gene_mask] = Wtot2

    stdG_a = torch.sqrt(EG2_a)
    stdG_b = torch.sqrt(EG2_b)
    stdG_a[stdG_a == 0] = 1
    stdG_b[stdG_b == 0] = 1

    # Compute gene-pair Z-scores
    if isinstance(C_gp, tuple):
        C_gp_0, C_gp_1 = C_gp
        z_0 = C_gp_0 / stdG_a
        z_1 = C_gp_1 / stdG_b
        mask = torch.abs(z_0) < torch.abs(z_1)
        Z_gp = torch.where(mask, z_0, z_1)
        EG2_gp = torch.where(mask, EG2_a, EG2_b)
    else:
        C_gp = C_gp
        z_a = C_gp / stdG_a
        z_b = C_gp / stdG_b
        mask = torch.abs(z_a) < torch.abs(z_b)
        Z_gp = torch.where(mask, z_a, z_b)
        EG2_gp = torch.where(mask, EG2_a, EG2_b)

    # Compute metabolite-level expected variance
    EG2_m = compute_metabolite_cs(EG2_gp, gene_pair_dict, interacting_cell_scores=False)
    if not isinstance(EG2_m, torch.Tensor):
        EG2_m = torch.tensor(EG2_m, device=device, dtype=torch.float64)

    stdG_m = torch.sqrt(EG2_m)
    stdG_m[stdG_m == 0] = 1

    # Compute metabolite Z-scores
    if isinstance(C_m, tuple):
        C_m_0, C_m_1 = C_m
        z_0 = C_m_0 / stdG_m
        z_1 = C_m_1 / stdG_m
        Z_m = torch.where(torch.abs(z_0) < torch.abs(z_1), z_0, z_1)
    else:
        Z_m = C_m / stdG_m

    return Z_gp, Z_m


def compute_p_results_old(
    ct_pair,
    cell_type_pairs,
    cs_gp,
    cs_m,
    gene_pairs_ind,
    gene_pairs_ind_per_ct_pair,
    Wtot2,
    eg2s_gp,
    cell_type_key,
    gene_pair_dict,
):
    i = cell_type_pairs.index(ct_pair)
    gene_pair_cor_ct = (
        (cs_gp[0][i, :], cs_gp[1][i, :]) if isinstance(cs_gp, tuple) else cs_gp[i, :]
    )
    C_m = (cs_m[0][i, :], cs_m[1][i, :]) if isinstance(cs_m, tuple) else cs_m[i, :]
    gene_pairs_ind_ct_pair = gene_pairs_ind_per_ct_pair[
        ct_pair
    ]  # If we consider all the gene pairs (irrespective of the cell type pair) use 'gene_pairs_ind' directly

    eg2s_a, eg2s_b = eg2s_gp
    C_gp = []
    EG2_a = []
    EG2_b = []
    for gene_pair_ind_ct_pair in gene_pairs_ind_ct_pair:
        idx = gene_pairs_ind.index(gene_pair_ind_ct_pair)
        g1_ind, g2_ind = gene_pair_ind_ct_pair
        lc_gp = (
            (gene_pair_cor_ct[0][idx], gene_pair_cor_ct[1][idx])
            if isinstance(gene_pair_cor_ct, tuple)
            else gene_pair_cor_ct[idx]
        )
        eg2_a = eg2_b = Wtot2[i]
        # if g1_ind == g2_ind:
        #     eg2_a = eg2_b = Wtot2[i]
        # else:
        #     eg2_a = eg2s_a[i, idx]
        #     eg2_b = eg2s_b[i, idx]
        C_gp.append(lc_gp)
        EG2_a.append(eg2_a)
        EG2_b.append(eg2_b)

    # Gene pairs

    EG = [0 for i in range(len(gene_pairs_ind_ct_pair))]

    stdG_a = [(EG2_a[i] - EG[i] ** 2) ** 0.5 for i in range(len(gene_pairs_ind_ct_pair))]
    stdG_a = [1 if stdG_a[i] == 0 else stdG_a[i] for i in range(len(stdG_a))]

    stdG_b = [(EG2_b[i] - EG[i] ** 2) ** 0.5 for i in range(len(gene_pairs_ind_ct_pair))]
    stdG_b = [1 if stdG_b[i] == 0 else stdG_b[i] for i in range(len(stdG_b))]

    if isinstance(C_gp[0], tuple):
        Z_gp_a = [(C_gp[i][0] - EG[i]) / stdG_a[i] for i in range(len(gene_pairs_ind_ct_pair))]
        Z_gp_b = [(C_gp[i][1] - EG[i]) / stdG_b[i] for i in range(len(gene_pairs_ind_ct_pair))]
    else:
        Z_gp_a = [(C_gp[i] - EG[i]) / stdG_a[i] for i in range(len(gene_pairs_ind_ct_pair))]
        Z_gp_b = [(C_gp[i] - EG[i]) / stdG_b[i] for i in range(len(gene_pairs_ind_ct_pair))]

    EG2_m_a = compute_metabolite_cs_old(
        np.array(EG2_a),
        cell_type_key=None,
        gene_pair_dict=gene_pair_dict,
        interacting_cell_scores=False,
    )
    EG2_m_b = compute_metabolite_cs_old(
        np.array(EG2_b),
        cell_type_key=None,
        gene_pair_dict=gene_pair_dict,
        interacting_cell_scores=False,
    )

    # Metabolites

    EG = [0 for i in range(len(gene_pair_dict.keys()))]

    stdG_a = [(EG2_m_a[i] - EG[i] ** 2) ** 0.5 for i in range(len(gene_pair_dict.keys()))]
    stdG_a = [1 if stdG_a[i] == 0 else stdG_a[i] for i in range(len(stdG_a))]

    stdG_b = [(EG2_m_b[i] - EG[i] ** 2) ** 0.5 for i in range(len(gene_pair_dict.keys()))]
    stdG_b = [1 if stdG_b[i] == 0 else stdG_b[i] for i in range(len(stdG_b))]

    if isinstance(C_m, tuple):
        Z_m_a = [(C_m[0][i] - EG[i]) / stdG_a[i] for i in range(len(gene_pair_dict.keys()))]
        Z_m_b = [(C_m[1][i] - EG[i]) / stdG_b[i] for i in range(len(gene_pair_dict.keys()))]
    else:
        Z_m_a = [(C_m[i] - EG[i]) / stdG_a[i] for i in range(len(gene_pair_dict.keys()))]
        Z_m_b = [(C_m[i] - EG[i]) / stdG_b[i] for i in range(len(gene_pair_dict.keys()))]

    return (C_gp, Z_gp_a, Z_gp_b, C_m, Z_m_a, Z_m_b)


def compute_p_int_cell_results(
    ct_pair,
    cell_type_pairs,
    cs_gp,
    cs_m,
    gene_pairs_ind,
    gene_pairs_ind_per_ct_pair,
    Wtot2,
    eg2s_gp,
    cell_type_key,
    gene_pair_dict,
):
    i = cell_type_pairs.index(ct_pair)
    gene_pair_cor_ct = cs_gp[i, :, :]
    C_m = cs_m[i, :, :]
    gene_pairs_ind_ct_pair = gene_pairs_ind_per_ct_pair[
        ct_pair
    ]  # If we consider all the gene pairs (irrespective of the cell type pair) use 'gene_pairs_ind' directly

    eg2s_a, eg2s_b = eg2s_gp
    C_gp = []
    EG2_a = []
    EG2_b = []
    for gene_pair_ind_ct_pair in gene_pairs_ind_ct_pair:
        idx = gene_pairs_ind.index(gene_pair_ind_ct_pair)
        g1_ind, g2_ind = gene_pair_ind_ct_pair
        lc_gp = gene_pair_cor_ct[:, idx]
        if g1_ind == g2_ind:
            eg2_a = eg2_b = Wtot2[i, :]
        else:
            eg2_a = (
                eg2s_a[i, :, g1_ind]
                if type(g1_ind) is not list
                else np.max(eg2s_a[i, :, g1_ind], axis=0)
            )
            eg2_b = (
                eg2s_b[i, :, g2_ind]
                if type(g2_ind) is not list
                else np.max(eg2s_b[i, :, g2_ind], axis=0)
            )
        C_gp.append(lc_gp)
        EG2_a.append(eg2_a)
        EG2_b.append(eg2_b)
    C_gp = np.column_stack(C_gp)
    EG2_a = np.column_stack(EG2_a)
    EG2_b = np.column_stack(EG2_b)

    # Gene pairs

    EG = np.zeros(C_gp.shape)

    stdG_a = (EG2_a - EG**2) ** 0.5
    stdG_a[stdG_a == 0] = 1

    stdG_b = (EG2_b - EG**2) ** 0.5
    stdG_b[stdG_b == 0] = 1

    Z_gp = np.where(
        np.abs((C_gp - EG) / stdG_a) < np.abs((C_gp - EG) / stdG_b),
        (C_gp - EG) / stdG_a,
        (C_gp - EG) / stdG_b,
    )

    EG2_gp = np.where(np.abs((C_gp - EG) / stdG_a) < np.abs((C_gp - EG) / stdG_b), EG2_a, EG2_b)

    EG2_m = compute_metabolite_cs(
        EG2_gp, cell_type_key=None, gene_pair_dict=gene_pair_dict, interacting_cell_scores=True
    )

    # Metabolites

    EG = np.zeros(C_m.shape)

    stdG = (EG2_m - EG**2) ** 0.5
    stdG[stdG == 0] = 1

    Z_m = (C_m - EG) / stdG

    return (C_gp, Z_gp, C_m, Z_m)


def compute_p_int_cell_results_no_ct(C_gp, C_m, gene_pairs_ind, Wtot2, eg2s_gp, gene_pair_dict):

    device = Wtot2.device

    # Convert indices
    same_gene_mask = torch.tensor(
        [
            (isinstance(g1, int) and isinstance(g2, int) and g1 == g2)
            or (isinstance(g1, list) and isinstance(g2, list) and sorted(g1) == sorted(g2))
            for g1, g2 in gene_pairs_ind
        ],
        device=device,
    )

    # Unpack second moments
    EG2_a = eg2s_gp[0].clone()
    EG2_b = eg2s_gp[1].clone()
    EG2_a[:, same_gene_mask] = Wtot2
    EG2_b[:, same_gene_mask] = Wtot2

    stdG_a = torch.sqrt(EG2_a)
    stdG_b = torch.sqrt(EG2_b)
    stdG_a[stdG_a == 0] = 1
    stdG_b[stdG_b == 0] = 1

    # Compute gene-pair Z-scores
    if isinstance(C_gp, tuple):
        C_gp_0, C_gp_1 = C_gp
        z_0 = C_gp_0 / stdG_a
        z_1 = C_gp_1 / stdG_b
        mask = torch.abs(z_0) < torch.abs(z_1)
        Z_gp = torch.where(mask, z_0, z_1)
        EG2_gp = torch.where(mask, EG2_a, EG2_b)
    else:
        C_gp = C_gp
        z_a = C_gp / stdG_a
        z_b = C_gp / stdG_b
        mask = torch.abs(z_a) < torch.abs(z_b)
        Z_gp = torch.where(mask, z_a, z_b)
        EG2_gp = torch.where(mask, EG2_a, EG2_b)

    # Compute metabolite-level expected variance
    EG2_m = compute_metabolite_cs(EG2_gp, gene_pair_dict, interacting_cell_scores=True)
    if not isinstance(EG2_m, torch.Tensor):
        EG2_m = torch.tensor(EG2_m, device=device, dtype=torch.float64)

    stdG_m = torch.sqrt(EG2_m)
    stdG_m[stdG_m == 0] = 1

    # Compute metabolite Z-scores
    if isinstance(C_m, tuple):
        C_m_0, C_m_1 = C_m
        z_0 = C_m_0 / stdG_m
        z_1 = C_m_1 / stdG_m
        Z_m = torch.where(torch.abs(z_0) < torch.abs(z_1), z_0, z_1)
    else:
        Z_m = C_m / stdG_m

    return (Z_gp, Z_m)


def compute_np_results(
    ct_pair,
    cell_type_pairs,
    cs_gp,
    cs_m,
    pvals_gp,
    pvals_m,
    gene_pair_dict,
    gene_pairs_ind,
    gene_pairs_ind_per_ct_pair,
):
    i = cell_type_pairs.index(ct_pair)
    gene_pair_cor_gp_ct = cs_gp[i, :]
    C_m = cs_m[i, :]
    pvals_gp_a, pvals_gp_b = pvals_gp
    pvals_gp_a_ct = pvals_gp_a[i, :]
    pvals_gp_b_ct = pvals_gp_b[i, :]
    pvals_m_a, pvals_m_b = pvals_m
    p_values_m_a = pvals_m_a[i, :]
    gene_pairs_ind_ct_pair = gene_pairs_ind_per_ct_pair[
        ct_pair
    ]  # If we consider all the gene pairs (irrespective of the cell type pair) use 'gene_pairs_ind' directly

    C_gp = []
    p_values_gp_a = []
    p_values_gp_b = []
    for gene_pair_ind_ct_pair in gene_pairs_ind_ct_pair:
        idx = gene_pairs_ind.index(gene_pair_ind_ct_pair)
        lc_gp = gene_pair_cor_gp_ct[idx]
        p_value_gp_a = pvals_gp_a_ct[idx]
        p_value_gp_b = pvals_gp_b_ct[idx]

        C_gp.append(lc_gp.reshape(1))
        p_values_gp_a.append(p_value_gp_a.reshape(1))
        p_values_gp_b.append(p_value_gp_b.reshape(1))

    C_gp = list(np.concatenate(C_gp))
    p_values_gp_a = list(np.concatenate(p_values_gp_a))
    p_values_gp_b = list(np.concatenate(p_values_gp_b))

    # C_m = compute_metabolite_cs(np.array(C_gp), cell_type_key=None, gene_pair_dict=gene_pair_dict, interacting_cell_scores=False)

    return (C_gp, p_values_gp_a, p_values_gp_b, C_m, p_values_m_a, p_values_m_a)


def get_ct_cell_communication_results(
    adata,
    genes,
    cells,
    layer_key_p_test,
    layer_key_np_test,
    model,
    cell_types,
    cell_type_pairs,
    D,
    test,
    device,
):

    gene_pairs_ind_per_ct_pair = adata.uns["gene_pairs_ind_per_ct_pair"]
    gene_pair_dict = adata.uns["gene_pair_dict"]
    genes = adata.uns["genes"]

    sample_specific = "sample_key" in adata.uns

    if isinstance(D, np.ndarray):
        D = torch.tensor(D, dtype=torch.float64, device=device)

    def idx_to_gene(idx):
        return [genes[i] for i in idx] if isinstance(idx, list) else genes[idx]

    records = [
        {
            "Cell Type 1": ct1,
            "Cell Type 2": ct2,
            "Gene 1": idx_to_gene(gp[0]),
            "Gene 2": idx_to_gene(gp[1]),
        }
        for (ct1, ct2), gp_list in gene_pairs_ind_per_ct_pair.items()
        for gp in gp_list
    ]
    cell_com_df_gp = pd.DataFrame.from_records(records)

    # Generate metabolite interaction table
    ct_pairs = list(gene_pairs_ind_per_ct_pair.keys())
    metabolites = list(gene_pair_dict.keys())
    cell_com_df_m = pd.DataFrame(
        [
            {"Cell Type 1": ct1, "Cell Type 2": ct2, "metabolite": m}
            for (ct1, ct2), m in itertools.product(ct_pairs, metabolites)
        ]
    )

    if test in ["parametric", "both"]:
        suffix = "p"
        # Gene pair
        c_values = adata.uns["ct_ccc_results"][suffix]["gp"]["cs"]
        z_values = adata.uns["ct_ccc_results"][suffix]["gp"]["Z"]
        p_values = adata.uns["ct_ccc_results"][suffix]["gp"]["Z_pval"]
        fdr_values = adata.uns["ct_ccc_results"][suffix]["gp"]["Z_FDR"]
        cell_com_df_gp[f"C_{suffix}"] = c_values.flatten()
        cell_com_df_gp["Z"] = z_values.flatten()
        cell_com_df_gp["Z_pval"] = p_values.flatten()
        cell_com_df_gp["Z_FDR"] = fdr_values.flatten()

        counts = counts_from_anndata(adata[:, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)
        counts_std = standardize_ct_counts(
            adata, counts, model, num_umi, sample_specific, cell_types
        )

        c_values_norm = normalize_ct_values(
            counts_std, cell_types, cell_type_pairs, gene_pairs_ind_per_ct_pair, c_values, D
        )
        adata.uns["ct_ccc_results"][suffix]["gp"]["cs_norm"] = c_values_norm.cpu().numpy()
        cell_com_df_gp[f"C_norm_{suffix}"] = c_values_norm.cpu().numpy().flatten()

        # Metabolite
        c_values = adata.uns["ct_ccc_results"][suffix]["m"]["cs"]
        z_values = adata.uns["ct_ccc_results"][suffix]["m"]["Z"]
        p_values = adata.uns["ct_ccc_results"][suffix]["m"]["Z_pval"]
        fdr_values = adata.uns["ct_ccc_results"][suffix]["m"]["Z_FDR"]
        cell_com_df_m[f"C_{suffix}"] = c_values.flatten()
        cell_com_df_m["Z"] = z_values.flatten()
        cell_com_df_m["Z_pval"] = p_values.flatten()
        cell_com_df_m["Z_FDR"] = fdr_values.flatten()

    if test in ["non-parametric", "both"]:
        suffix = "np"
        # Gene pair
        c_values = adata.uns["ct_ccc_results"][suffix]["gp"]["cs"]
        p_values = adata.uns["ct_ccc_results"][suffix]["gp"]["pval"]
        fdr_values = adata.uns["ct_ccc_results"][suffix]["gp"]["FDR"]
        cell_com_df_gp[f"C_{suffix}"] = c_values.flatten()
        cell_com_df_gp[f"pval_{suffix}"] = p_values.flatten()
        cell_com_df_gp[f"FDR_{suffix}"] = fdr_values.flatten()

        counts = counts_from_anndata(adata[:, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        if adata.uns.get("center_counts_for_np_test", False):
            num_umi = counts.sum(dim=0)
            counts = standardize_ct_counts(
                adata, counts, model, num_umi, sample_specific, cell_types
            )

        c_values_norm = normalize_ct_values(
            counts, cell_types, cell_type_pairs, gene_pairs_ind_per_ct_pair, c_values, D
        )
        adata.uns["ct_ccc_results"][suffix]["gp"]["cs_norm"] = c_values_norm.cpu().numpy()
        cell_com_df_gp[f"C_norm_{suffix}"] = c_values_norm.cpu().numpy().flatten()

        # Metabolite
        c_values = adata.uns["ct_ccc_results"][suffix]["m"]["cs"]
        p_values = adata.uns["ct_ccc_results"][suffix]["m"]["pval"]
        fdr_values = adata.uns["ct_ccc_results"][suffix]["m"]["FDR"]
        cell_com_df_m[f"C_{suffix}"] = c_values.flatten()
        cell_com_df_m[f"pval_{suffix}"] = p_values.flatten()
        cell_com_df_m[f"FDR_{suffix}"] = fdr_values.flatten()

    adata.uns["ct_ccc_results"]["cell_com_df_gp"] = cell_com_df_gp
    adata.uns["ct_ccc_results"]["cell_com_df_m"] = cell_com_df_m

    return


def get_cell_communication_results_old(
    adata,
    genes,
    cells,
    layer_key_p_test,
    layer_key_np_test,
    model,
    cell_types,
    cell_type_pairs,
    D,
    test,
):

    gene_pairs_ind_per_ct_pair = adata.uns["gene_pairs_ind_per_ct_pair"]
    gene_pair_dict = adata.uns["gene_pair_dict"]
    genes = adata.uns["genes"]

    def map_to_genes(x):
        if isinstance(x, list):
            return [genes[i] for i in x]
        else:
            return genes[x]

    cell_com_df_gp = (
        pd.DataFrame.from_dict(gene_pairs_ind_per_ct_pair, orient="index")
        .stack()
        .to_frame()
        .reset_index()
    )
    cell_com_df_gp = cell_com_df_gp.drop(["level_1"], axis=1)
    cell_com_df_gp = cell_com_df_gp.rename(columns={"level_0": "cell_type_pair", 0: "gene_pair"})
    cell_com_df_gp["Cell Type 1"], cell_com_df_gp["Cell Type 2"] = zip(
        *cell_com_df_gp["cell_type_pair"], strict=False
    )
    cell_com_df_gp["Gene 1"], cell_com_df_gp["Gene 2"] = zip(
        *cell_com_df_gp["gene_pair"], strict=False
    )
    cell_com_df_gp["Gene 1"] = cell_com_df_gp["Gene 1"].apply(map_to_genes)
    cell_com_df_gp["Gene 2"] = cell_com_df_gp["Gene 2"].apply(map_to_genes)
    cell_com_df_gp = cell_com_df_gp.drop(["cell_type_pair", "gene_pair"], axis=1)

    ct_pair_metab = list(
        itertools.product(gene_pairs_ind_per_ct_pair.keys(), gene_pair_dict.keys())
    )
    cell_com_df_m = pd.DataFrame(ct_pair_metab, columns=["cell_type_pair", "metabolite"])
    cell_com_df_m["Cell Type 1"], cell_com_df_m["Cell Type 2"] = zip(
        *cell_com_df_m["cell_type_pair"], strict=False
    )
    cell_com_df_m = cell_com_df_m.drop(["cell_type_pair"], axis=1)

    if test in ["parametric", "both"]:
        # Gene pair
        c_values = adata.uns["ct_ccc_results"]["p"]["gp"]["cs"]
        cell_com_df_gp["C_p"] = c_values.flatten()
        z_values_a = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_a"]
        cell_com_df_gp["Z_a"] = z_values_a.flatten()
        z_values_b = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_b"]
        cell_com_df_gp["Z_b"] = z_values_b.flatten()
        p_values_a = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_pval_a"]
        cell_com_df_gp["Z_pval_a"] = p_values_a.flatten()
        p_values_b = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_pval_b"]
        cell_com_df_gp["Z_pval_b"] = p_values_b.flatten()
        FDR_values_a = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_FDR_a"]
        cell_com_df_gp["Z_FDR_a"] = FDR_values_a.flatten()
        FDR_values_b = adata.uns["ct_ccc_results"]["p"]["gp"]["Z_FDR_b"]
        cell_com_df_gp["Z_FDR_b"] = FDR_values_b.flatten()

        counts = counts_from_anndata(adata[cells, genes], layer_key_p_test, dense=True)
        num_umi = counts.sum(axis=0)
        counts_std = counts_std = create_centered_counts_ct(counts, model, num_umi, cell_types)
        counts_std = np.nan_to_num(counts_std)

        c_values_norm = normalize_values_old(
            counts_std, cell_types, cell_type_pairs, gene_pairs_ind_per_ct_pair, c_values, D
        )
        adata.uns["ct_ccc_results"]["p"]["gp"]["cs_norm"] = c_values_norm
        cell_com_df_gp["C_norm_p"] = c_values_norm.flatten()

        # Metabolite
        c_values = adata.uns["ct_ccc_results"]["p"]["m"]["cs"]
        cell_com_df_m["C_p"] = c_values.flatten()
        z_values_a = adata.uns["ct_ccc_results"]["p"]["m"]["Z_a"]
        cell_com_df_m["Z_a"] = z_values_a.flatten()
        z_values_b = adata.uns["ct_ccc_results"]["p"]["m"]["Z_b"]
        cell_com_df_m["Z_b"] = z_values_b.flatten()
        p_values_a = adata.uns["ct_ccc_results"]["p"]["m"]["Z_pval_a"]
        cell_com_df_m["Z_pval_a"] = p_values_a.flatten()
        p_values_b = adata.uns["ct_ccc_results"]["p"]["m"]["Z_pval_b"]
        cell_com_df_m["Z_pval_b"] = p_values_b.flatten()
        FDR_values_a = adata.uns["ct_ccc_results"]["p"]["m"]["Z_FDR_a"]
        cell_com_df_m["Z_FDR_a"] = FDR_values_a.flatten()
        FDR_values_b = adata.uns["ct_ccc_results"]["p"]["m"]["Z_FDR_b"]
        cell_com_df_m["Z_FDR_b"] = FDR_values_b.flatten()

    if test in ["non-parametric", "both"]:
        # Gene pair
        c_values = adata.uns["ct_ccc_results"]["np"]["gp"]["cs"]
        cell_com_df_gp["C_np"] = c_values.flatten()
        p_values_a = adata.uns["ct_ccc_results"]["np"]["gp"]["pval_a"]
        cell_com_df_gp["pval_np_a"] = p_values_a.flatten()
        p_values_b = adata.uns["ct_ccc_results"]["np"]["gp"]["pval_b"]
        cell_com_df_gp["pval_np_b"] = p_values_b.flatten()
        FDR_values_a = adata.uns["ct_ccc_results"]["np"]["gp"]["FDR_a"]
        cell_com_df_gp["FDR_np_a"] = FDR_values_a.flatten()
        FDR_values_b = adata.uns["ct_ccc_results"]["np"]["gp"]["FDR_b"]
        cell_com_df_gp["FDR_np_b"] = FDR_values_b.flatten()

        counts = counts_from_anndata(adata[:, genes], layer_key_np_test, dense=True)
        if adata.uns["center_counts_for_np_test"]:
            num_umi = counts.sum(axis=0)
            counts = create_centered_counts(counts, model, num_umi)
            counts = np.nan_to_num(counts)
        c_values_norm = normalize_values_old(
            counts, cell_types, cell_type_pairs, gene_pairs_ind_per_ct_pair, c_values, D
        )
        adata.uns["ct_ccc_results"]["np"]["gp"]["cs_norm"] = c_values_norm
        cell_com_df_gp["C_norm_np"] = c_values_norm.flatten()

        # Metabolite
        c_values = adata.uns["ct_ccc_results"]["np"]["m"]["cs"]
        cell_com_df_m["C_np"] = c_values.flatten()
        p_values_a = adata.uns["ct_ccc_results"]["np"]["m"]["pval_a"]
        cell_com_df_m["pval_np_a"] = p_values_a.flatten()
        p_values_b = adata.uns["ct_ccc_results"]["np"]["m"]["pval_b"]
        cell_com_df_m["pval_np_b"] = p_values_b.flatten()
        FDR_values_a = adata.uns["ct_ccc_results"]["np"]["m"]["FDR_a"]
        cell_com_df_m["FDR_np_a"] = FDR_values_a.flatten()
        FDR_values_b = adata.uns["ct_ccc_results"]["np"]["m"]["FDR_b"]
        cell_com_df_m["FDR_np_b"] = FDR_values_b.flatten()

    adata.uns["ct_ccc_results"]["cell_com_df_gp"] = cell_com_df_gp
    adata.uns["ct_ccc_results"]["cell_com_df_m"] = cell_com_df_m

    return


def get_cell_communication_results(
    adata,
    genes,
    layer_key_p_test,
    layer_key_np_test,
    model,
    D,
    test,
    device,
):

    gene_pairs = adata.uns["gene_pairs"]
    gene_pairs_ind = adata.uns["gene_pairs_ind"]
    gene_pair_dict = adata.uns["gene_pair_dict"]

    sample_specific = "sample_key" in adata.uns

    if isinstance(D, np.ndarray):
        D = torch.tensor(D, dtype=torch.float64, device=device)

    # Initialize dataframes
    cell_com_df_gp = pd.DataFrame(gene_pairs, columns=["Gene 1", "Gene 2"])
    cell_com_df_m = pd.DataFrame({"Metabolite": list(gene_pair_dict.keys())})

    if test in ["parametric", "both"]:
        suffix = "p"
        # Gene pair
        c_values = adata.uns["ccc_results"][suffix]["gp"]["cs"]
        z_values = adata.uns["ccc_results"][suffix]["gp"]["Z"]
        p_values = adata.uns["ccc_results"][suffix]["gp"]["Z_pval"]
        fdr_values = adata.uns["ccc_results"][suffix]["gp"]["Z_FDR"]
        cell_com_df_gp[f"C_{suffix}"] = c_values
        cell_com_df_gp["Z"] = z_values
        cell_com_df_gp["Z_pval"] = p_values
        cell_com_df_gp["Z_FDR"] = fdr_values

        counts = counts_from_anndata(adata[:, genes], layer_key_p_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        num_umi = counts.sum(dim=0)
        _lazy_import_hotspot()
        counts_std = standardize_counts(adata, counts, model, num_umi, sample_specific)

        c_values_norm = normalize_values(counts_std, gene_pairs_ind, c_values, D)
        adata.uns["ccc_results"][suffix]["gp"]["cs_norm"] = c_values_norm.cpu().numpy()
        cell_com_df_gp[f"C_norm_{suffix}"] = c_values_norm.cpu().numpy()

        # Metabolite
        c_values = adata.uns["ccc_results"][suffix]["m"]["cs"]
        z_values = adata.uns["ccc_results"][suffix]["m"]["Z"]
        p_values = adata.uns["ccc_results"][suffix]["m"]["Z_pval"]
        fdr_values = adata.uns["ccc_results"][suffix]["m"]["Z_FDR"]
        cell_com_df_m[f"C_{suffix}"] = c_values
        cell_com_df_m["Z"] = z_values
        cell_com_df_m["Z_pval"] = p_values
        cell_com_df_m["Z_FDR"] = fdr_values

    if test in ["non-parametric", "both"]:
        suffix = "np"
        # Gene pair
        c_values = adata.uns["ccc_results"][suffix]["gp"]["cs"]
        p_values = adata.uns["ccc_results"][suffix]["gp"]["pval"]
        fdr_values = adata.uns["ccc_results"][suffix]["gp"]["FDR"]
        cell_com_df_gp[f"C_{suffix}"] = c_values
        cell_com_df_gp[f"pval_{suffix}"] = p_values
        cell_com_df_gp[f"FDR_{suffix}"] = fdr_values

        counts = counts_from_anndata(adata[:, genes], layer_key_np_test, dense=True)
        counts = torch.tensor(counts, dtype=torch.float64, device=device)
        if adata.uns.get("center_counts_for_np_test", False):
            num_umi = counts.sum(dim=0)
            _lazy_import_hotspot()
            counts = standardize_counts(adata, counts, model, num_umi, sample_specific)

        c_values_norm = normalize_values(counts, gene_pairs_ind, c_values, D)
        adata.uns["ccc_results"][suffix]["gp"]["cs_norm"] = c_values_norm.cpu().numpy()
        cell_com_df_gp[f"C_norm_{suffix}"] = c_values_norm.cpu().numpy()

        # Metabolite
        c_values = adata.uns["ccc_results"][suffix]["m"]["cs"]
        p_values = adata.uns["ccc_results"][suffix]["m"]["pval"]
        fdr_values = adata.uns["ccc_results"][suffix]["m"]["FDR"]
        cell_com_df_m[f"C_{suffix}"] = c_values
        cell_com_df_m[f"pval_{suffix}"] = p_values
        cell_com_df_m[f"FDR_{suffix}"] = fdr_values

    adata.uns["ccc_results"]["cell_com_df_gp"] = cell_com_df_gp
    adata.uns["ccc_results"]["cell_com_df_m"] = cell_com_df_m

    return


def normalize_ct_values(
    counts,
    cell_types,
    cell_type_pairs,
    gene_pairs_per_ct_pair_ind,
    lcs,
    D,
):

    if isinstance(cell_types, pd.Series):
        cell_types = cell_types.values

    if isinstance(lcs, np.ndarray):
        lcs = torch.tensor(lcs, dtype=counts.dtype, device=counts.device)

    c_values_norm = torch.empty_like(lcs, dtype=counts.dtype, device=counts.device)

    for i, ct_pair in enumerate(cell_type_pairs):
        ct_t, _ = ct_pair

        ct_mask = cell_types == ct_t
        if isinstance(ct_mask, np.ndarray):
            ct_mask = torch.tensor(ct_mask, device=counts.device)

        counts_ct = counts[:, ct_mask]
        D_ct = D[i][ct_mask]
        gene_pairs_ind = gene_pairs_per_ct_pair_ind[ct_pair]

        lc_maxs = compute_max_cs(D_ct, counts_ct, gene_pairs_ind)
        lc_maxs = torch.where(lc_maxs == 0, torch.tensor(1.0, device=counts.device), lc_maxs)

        c_values = lcs[i] if lcs.ndim == 2 else lcs[i : i + 1]  # allow 1D or 2D lcs
        c_values_norm[i] = c_values / lc_maxs
        c_values_norm[i] = torch.where(
            torch.isinf(c_values_norm[i]),
            torch.tensor(1.0, device=counts.device),
            c_values_norm[i],
        )

    return c_values_norm


def normalize_values(counts, gene_pairs_ind, lcs, D):
    """Normalize communication scores (lcs) using maximum possible score estimates."""
    lc_maxs = compute_max_cs(D, counts, gene_pairs_ind)
    lc_maxs = torch.where(lc_maxs == 0, torch.tensor(1.0, device=lc_maxs.device), lc_maxs)
    if isinstance(lcs, np.ndarray):
        lcs = torch.tensor(lcs, dtype=lc_maxs.dtype, device=lc_maxs.device)
    c_values_norm = lcs / lc_maxs
    c_values_norm = torch.where(
        torch.isinf(c_values_norm), torch.tensor(1.0, device=c_values_norm.device), c_values_norm
    )
    return c_values_norm


def normalize_values_old(
    counts,
    cell_types,
    cell_type_pairs,
    gene_pairs_per_ct_pair_ind,
    lcs,
    D,
):

    c_values_norm = np.zeros(lcs.shape)
    for i in range(len(cell_type_pairs)):
        ct_pair = cell_type_pairs[i]
        ct_t, ct_u = ct_pair
        cell_type_t_mask = [ct == ct_t for ct in cell_types]
        counts_ct_t = counts[:, cell_type_t_mask]
        D_ct_t = D[i][cell_type_t_mask]
        gene_pairs_ind = gene_pairs_per_ct_pair_ind[ct_pair]
        lc_maxs = compute_max_cs_old(D_ct_t, counts_ct_t, gene_pairs_ind)
        lc_maxs[lc_maxs == 0] = 1
        c_values_norm[i] = lcs[i] / lc_maxs
        c_values_norm[i][np.isinf(c_values_norm[i])] = 1

    return c_values_norm


def compute_max_cs(node_degrees, counts, gene_pairs_ind):
    """Compute max communication scores per gene pair."""
    result = torch.empty(len(gene_pairs_ind), dtype=counts.dtype, device=counts.device)

    for i, (g1, _) in enumerate(gene_pairs_ind):
        if isinstance(g1, list):
            vals = counts[g1].mean(dim=0)
        else:
            vals = counts[g1]
        result[i] = compute_max_cs_gp(vals, node_degrees)

    return result


def compute_max_cs_old(node_degrees, counts, gene_pairs_ind):

    result = np.zeros(len(gene_pairs_ind))
    for i, gene_pair_ind in enumerate(gene_pairs_ind):
        vals = (
            counts[gene_pair_ind[0]]
            if type(gene_pair_ind[0]) is not list
            else np.mean(counts[gene_pair_ind[0]], axis=0)
        )
        result[i] = compute_max_cs_gp_old(vals, node_degrees)

    return result


def compute_max_cs_gp(vals, node_degrees):
    """Compute max communication score for a single gene (vector)."""
    return 0.5 * torch.sum(node_degrees * vals**2)


@jit(nopython=True)
def compute_max_cs_gp_old(vals, node_degrees):
    tot = 0.0

    for i in range(node_degrees.size):
        tot += node_degrees[i] * (vals[i] ** 2)

    return tot / 2


def center_ct_counts_torch(counts, num_umi, model, cell_types):
    """
    counts: Tensor [genes, cells]
    num_umi: Tensor [cells]
    model: 'bernoulli', 'danb', 'normal', or 'none'

    Returns
    -------
        Centered counts within cell types: Tensor [genes, cells]
    """
    # Binarize if using Bernoulli
    if model == "bernoulli":
        counts = (counts > 0).double()
        mu, var, _ = models.apply_model_per_cell_type(
            models.bernoulli_model_torch, counts, num_umi, cell_types
        )
    elif model == "danb":
        mu, var, _ = models.apply_model_per_cell_type(
            models.danb_model_torch, counts, num_umi, cell_types
        )
    elif model == "normal":
        mu, var, _ = models.apply_model_per_cell_type(
            models.normal_model_torch, counts, num_umi, cell_types
        )
    elif model == "none":
        mu, var, _ = models.apply_model_per_cell_type(
            models.none_model_torch, counts, num_umi, cell_types
        )
    else:
        raise ValueError(f"Unsupported model type: {model}")

    # Avoid division by zero
    std = torch.sqrt(var)
    std[std == 0] = 1.0

    centered = (counts - mu) / std
    centered[centered == 0] = 0  # Optional: to match old behavior

    return centered


def create_centered_counts(counts, model, num_umi):
    """
    Creates a matrix of centered/standardized counts given
    the selected statistical model
    """
    out = np.zeros_like(counts, dtype="double")

    for i in range(out.shape[0]):
        vals_x = counts[i]

        out_x = create_centered_counts_row(vals_x, model, num_umi)

        out[i] = out_x

    return out


def create_centered_counts_ct(counts, model, num_umi, cell_types):
    """
    Creates a matrix of centered/standardized counts given
    the selected statistical model
    """
    out = np.zeros_like(counts, dtype="double")

    for i in range(out.shape[0]):
        vals_x = counts[i]

        out_x = create_centered_counts_row_ct(vals_x, model, num_umi, cell_types)

        out[i] = out_x

    return out


def create_centered_counts_row_ct(vals_x, model, num_umi, cell_types):
    if model == "bernoulli":
        vals_x = (vals_x > 0).astype("double")
        mu_x, var_x, x2_x = models.bernoulli_model(vals_x, num_umi)

    elif model == "danb":
        mu_x, var_x, x2_x = models.ct_danb_model(vals_x, num_umi, cell_types)

    elif model == "normal":
        mu_x, var_x, x2_x = models.normal_model(vals_x, num_umi)

    elif model == "none":
        mu_x, var_x, x2_x = models.none_model(vals_x, num_umi)

    else:
        raise Exception(f"Invalid Model: {model}")

    var_x[var_x == 0] = 1
    out_x = (vals_x - mu_x) / (var_x**0.5)
    out_x[out_x == 0] = 0

    return out_x


def create_centered_counts_row(vals_x, model, num_umi):
    if model == "bernoulli":
        vals_x = (vals_x > 0).astype("double")
        mu_x, var_x, x2_x = models.bernoulli_model(vals_x, num_umi)

    elif model == "danb":
        mu_x, var_x, x2_x = models.danb_model(vals_x, num_umi)

    elif model == "normal":
        mu_x, var_x, x2_x = models.normal_model(vals_x, num_umi)

    elif model == "none":
        mu_x, var_x, x2_x = models.none_model(vals_x, num_umi)

    else:
        raise Exception(f"Invalid Model: {model}")

    var_x[var_x == 0] = 1
    out_x = (vals_x - mu_x) / (var_x**0.5)
    out_x[out_x == 0] = 0

    return out_x


def compute_Z_scores_cellcom_p(
    ct_pair, cell_type_pairs, gene_pair_cor, gene_pairs_per_ct_pair_ind, Wtot2, eg2s
):
    i = cell_type_pairs.index(ct_pair)
    gene_pair_cor_ct = gene_pair_cor[i, :, :]
    gene_pairs_ind = gene_pairs_per_ct_pair_ind[
        ct_pair
    ]  # If we consider all the gene pairs (irrespective of the cell type pair) use 'gene_pairs_ind' directly

    C = []
    EG2 = []
    for gene_pair_ind in gene_pairs_ind:
        g1_ind, g2_ind = gene_pair_ind
        lc = gene_pair_cor_ct[g1_ind, g2_ind]
        if g1_ind == g2_ind:
            eg2 = Wtot2[i]
        else:
            eg2 = eg2s[i][g1_ind]
        C.append(lc)
        EG2.append(eg2)

    EG = [0 for i in range(len(gene_pairs_ind))]

    stdG = [(EG2[i] - EG[i] ** 2) ** 0.5 for i in range(len(gene_pairs_ind))]
    stdG = [1 if stdG[i] == 0 else stdG[i] for i in range(len(stdG))]

    Z = [(C[i] - EG[i]) / stdG[i] for i in range(len(gene_pairs_ind))]

    return (C, Z)


def extract_results_cellcom_np(
    ct_pair, cell_type_pairs, gene_pair_cor, pvals, gene_pairs_per_ct_pair_ind
):
    i = cell_type_pairs.index(ct_pair)
    gene_pair_cor_ct = gene_pair_cor[i, :, :]
    pvals_ct = pvals[i, :, :]
    gene_pairs_ind = gene_pairs_per_ct_pair_ind[
        ct_pair
    ]  # If we consider all the gene pairs (irrespective of the cell type pair) use 'gene_pairs_ind' directly

    C = []
    p_values = []
    for gene_pair_ind in gene_pairs_ind:
        g1_ind, g2_ind = gene_pair_ind
        lc = gene_pair_cor_ct[g1_ind, g2_ind]
        p_value = pvals_ct[g1_ind, g2_ind]

        C.append(lc.reshape(1))
        p_values.append(p_value.reshape(1))

    C = list(np.concatenate(C))
    p_values = list(np.concatenate(p_values))

    return (C, p_values)


@njit
def expand_ct_pairs_cellcom(pairs, vals, N):
    out = [np.zeros((N, N)) for k in range(len(vals))]

    for k in range(len(out)):
        for i in range(len(pairs)):
            x = pairs[i, 0]
            y = pairs[i, 1]
            v = vals[k][i]

            out[k][x, y] = v

    return out


def compute_local_cov_pairs_max(node_degrees, counts):
    """
    For a Genes x Cells count matrix, compute the maximal pair-wise correlation
    between any two genes
    """
    N_GENES = counts.shape[0]

    gene_maxs = np.zeros(N_GENES)
    for i in range(N_GENES):
        gene_maxs[i] = compute_local_cov_max(counts[i].todense(), node_degrees)

    result = gene_maxs.reshape((-1, 1)) + gene_maxs.reshape((1, -1))
    result = result / 2
    return result


@jit(nopython=True)
def compute_local_cov_max(vals, node_degrees):
    tot = 0.0

    for i in range(node_degrees.size):
        tot += node_degrees[i] * (vals[i] ** 2)

    return tot / 2


def get_ct_pair_counts_and_weights(
    counts, weights, cell_type_pairs, cell_types, gene_pairs_per_ct_pair_ind
):

    c_nrow, c_ncol = counts.shape
    w_nrow, w_ncol = weights.shape
    n_ct_pairs = len(cell_type_pairs)

    extract_counts_weights_results = partial(
        extract_ct_pair_counts_weights,
        counts=counts,
        weights=weights,
        cell_type_pairs=cell_type_pairs,
        cell_types=cell_types,
        gene_pairs_per_ct_pair_ind=gene_pairs_per_ct_pair_ind,
    )
    results = list(map(extract_counts_weights_results, cell_type_pairs))

    c_new_data_t_all = [x[0] for x in results]
    c_new_coords_3d_t_all = [x[1] for x in results]
    c_new_coords_3d_t_all = np.hstack(c_new_coords_3d_t_all)
    c_new_data_t_all = np.concatenate(c_new_data_t_all)

    c_new_data_u_all = [x[2] for x in results]
    c_new_coords_3d_u_all = [x[3] for x in results]
    c_new_coords_3d_u_all = np.hstack(c_new_coords_3d_u_all)
    c_new_data_u_all = np.concatenate(c_new_data_u_all)

    w_new_data_all = [x[4] for x in results]
    w_new_coords_3d_all = [x[5] for x in results]
    w_new_coords_3d_all = np.hstack(w_new_coords_3d_all)
    w_new_data_all = np.concatenate(w_new_data_all)

    counts_ct_pairs_t = sparse.COO(
        c_new_coords_3d_t_all, c_new_data_t_all, shape=(n_ct_pairs, c_nrow, c_ncol)
    )
    counts_ct_pairs_u = sparse.COO(
        c_new_coords_3d_u_all, c_new_data_u_all, shape=(n_ct_pairs, c_nrow, c_ncol)
    )
    weights_ct_pairs = sparse.COO(
        w_new_coords_3d_all, w_new_data_all, shape=(n_ct_pairs, w_nrow, w_ncol)
    )

    return counts_ct_pairs_t, counts_ct_pairs_u, weights_ct_pairs


def get_ct_pair_counts_and_weights_null(
    counts_ct_pairs_t, counts_ct_pairs_u, weights, cell_type_pairs, cell_types, M
):

    n_cells = counts_ct_pairs_t.shape[2]
    cell_permutations = np.vstack([np.random.permutation(n_cells) for _ in range(M)])

    n_ct_pairs, c_nrow, c_ncol = counts_ct_pairs_t.shape
    w_nrow, w_ncol = weights.shape

    extract_counts_weights_results_null = partial(
        extract_ct_pair_counts_weights_null,
        permutations=cell_permutations,
        counts_ct_pairs_t=counts_ct_pairs_t,
        counts_ct_pairs_u=counts_ct_pairs_u,
        weights=weights,
        cell_type_pairs=cell_type_pairs,
        cell_types=cell_types,
    )
    results_null = list(map(extract_counts_weights_results_null, cell_permutations))

    c_null_data_t_all = [x[0] for x in results_null]
    c_null_coords_4d_t_all = [x[1] for x in results_null]
    c_null_coords_4d_t_all = np.hstack(c_null_coords_4d_t_all)
    c_null_data_t_all = np.concatenate(c_null_data_t_all)

    c_null_data_u_all = [x[2] for x in results_null]
    c_null_coords_4d_u_all = [x[3] for x in results_null]
    c_null_coords_4d_u_all = np.hstack(c_null_coords_4d_u_all)
    c_null_data_u_all = np.concatenate(c_null_data_u_all)

    w_null_data_all = [x[4] for x in results_null]
    w_null_coords_4d_all = [x[5] for x in results_null]
    w_null_coords_4d_all = np.hstack(w_null_coords_4d_all)
    w_null_data_all = np.concatenate(w_null_data_all)

    counts_ct_pairs_t_null = sparse.COO(
        c_null_coords_4d_t_all, c_null_data_t_all, shape=(M, n_ct_pairs, c_nrow, c_ncol)
    )
    counts_ct_pairs_u_null = sparse.COO(
        c_null_coords_4d_u_all, c_null_data_u_all, shape=(M, n_ct_pairs, c_nrow, c_ncol)
    )
    weights_ct_pairs_null = sparse.COO(
        w_null_coords_4d_all, w_null_data_all, shape=(M, n_ct_pairs, w_nrow, w_ncol)
    )

    return counts_ct_pairs_t_null, counts_ct_pairs_u_null, weights_ct_pairs_null


def extract_ct_pair_counts_weights(
    ct_pair, counts, weights, cell_type_pairs, cell_types, gene_pairs_per_ct_pair_ind
):

    i = cell_type_pairs.index(ct_pair)

    ct_t, ct_u = cell_type_pairs[i]
    gene_pairs_per_ct_pair_ind_i = gene_pairs_per_ct_pair_ind[(ct_t, ct_u)]
    ct_t_mask = cell_types.values == ct_t
    ct_t_mask_coords = np.argwhere(ct_t_mask)
    ct_u_mask = cell_types.values == ct_u
    ct_u_mask_coords = np.argwhere(ct_u_mask)

    ct_t_genes = np.unique([t[0] for t in gene_pairs_per_ct_pair_ind_i])
    ct_u_genes = np.unique([t[1] for t in gene_pairs_per_ct_pair_ind_i])

    c_old_coords = counts.coords
    w_old_coords = weights.coords

    # Counts

    c_row_coords_t, c_col_coords_t = np.meshgrid(ct_t_genes, ct_t_mask_coords, indexing="ij")
    c_row_coords_t = c_row_coords_t.ravel()
    c_col_coords_t = c_col_coords_t.ravel()
    c_new_coords_t = np.vstack((c_row_coords_t, c_col_coords_t))

    c_row_coords_u, c_col_coords_u = np.meshgrid(ct_u_genes, ct_u_mask_coords, indexing="ij")
    c_row_coords_u = c_row_coords_u.ravel()
    c_col_coords_u = c_col_coords_u.ravel()
    c_new_coords_u = np.vstack((c_row_coords_u, c_col_coords_u))

    c_matching_indices_t = np.where(np.all(np.isin(c_old_coords.T, c_new_coords_t.T), axis=1))[0]
    c_new_data_t = counts.data[c_matching_indices_t]
    c_new_coords_t = c_old_coords[:, c_matching_indices_t]

    c_matching_indices_u = np.where(np.all(np.isin(c_old_coords.T, c_new_coords_u.T), axis=1))[0]
    c_new_data_u = counts.data[c_matching_indices_u]
    c_new_coords_u = c_old_coords[:, c_matching_indices_u]

    c_coord_3d_t = np.full(c_new_coords_t.shape[1], fill_value=i)
    c_new_coords_3d_t = np.vstack((c_coord_3d_t, c_new_coords_t))

    c_coord_3d_u = np.full(c_new_coords_u.shape[1], fill_value=i)
    c_new_coords_3d_u = np.vstack((c_coord_3d_u, c_new_coords_u))

    # Weights

    w_row_coords, w_col_coords = np.meshgrid(ct_t_mask_coords, ct_u_mask_coords, indexing="ij")
    w_row_coords = w_row_coords.ravel()
    w_col_coords = w_col_coords.ravel()
    w_new_coords = np.vstack((w_row_coords, w_col_coords))

    w_matching_indices = np.where(np.all(np.isin(w_old_coords.T, w_new_coords.T), axis=1))[0]
    w_new_data = weights.data[w_matching_indices]
    w_new_coords = w_old_coords[:, w_matching_indices]

    w_coord_3d = np.full(w_new_coords.shape[1], fill_value=i)
    w_new_coords_3d = np.vstack((w_coord_3d, w_new_coords))

    return (
        c_new_data_t,
        c_new_coords_3d_t,
        c_new_data_u,
        c_new_coords_3d_u,
        w_new_data,
        w_new_coords_3d,
    )


def extract_ct_pair_counts_weights_null(
    permutation,
    permutations,
    counts_ct_pairs_t,
    counts_ct_pairs_u,
    weights,
    cell_type_pairs,
    cell_types,
):

    i = np.where(np.all(permutations == permutation, axis=1))[0]
    cell_types_perm = pd.Series(cell_types[permutation])

    counts_ct_pairs_t_perm = counts_ct_pairs_t[:, :, permutation]
    c_perm_coords_t = counts_ct_pairs_t_perm.coords
    c_perm_data_t = counts_ct_pairs_t_perm.data
    counts_ct_pairs_u_perm = counts_ct_pairs_u[:, :, permutation]
    c_perm_coords_u = counts_ct_pairs_u_perm.coords
    c_perm_data_u = counts_ct_pairs_u_perm.data

    extract_weights_results = partial(
        extract_ct_pair_weights,
        weights=weights,
        cell_type_pairs=cell_type_pairs,
        cell_types=cell_types_perm,
    )
    weights_results = list(map(extract_weights_results, cell_type_pairs))

    w_perm_data_all = [x[0] for x in weights_results]
    w_perm_coords_3d_all = [x[1] for x in weights_results]
    w_perm_coords_3d_all = np.hstack(w_perm_coords_3d_all)
    w_perm_data_all = np.concatenate(w_perm_data_all)

    c_coord_4d_t = np.full(c_perm_coords_t.shape[1], fill_value=i)
    c_perm_coords_4d_t = np.vstack((c_coord_4d_t, c_perm_coords_t))

    c_coord_4d_u = np.full(c_perm_coords_u.shape[1], fill_value=i)
    c_perm_coords_4d_u = np.vstack((c_coord_4d_u, c_perm_coords_u))

    w_coord_4d_u = np.full(w_perm_coords_3d_all.shape[1], fill_value=i)
    w_perm_coords_4d_u = np.vstack((w_coord_4d_u, w_perm_coords_3d_all))

    return (
        c_perm_data_t,
        c_perm_coords_4d_t,
        c_perm_data_u,
        c_perm_coords_4d_u,
        w_perm_data_all,
        w_perm_coords_4d_u,
    )


def compute_interaction_module_correlation(
    adata: AnnData,
    cor_method: Literal["pearson"] | Literal["spearman"] | None = "pearson",
    test: Literal["parametric"] | Literal["non-parametric"] | None = None,
    interaction_type: Literal["metabolite"] | Literal["gene_pair"] | None = "metabolite",
    only_sig_values: bool | None = False,
    normalize_values: bool | None = False,
    use_FDR: bool | None = True,
    use_super_modules: bool | None = False,
):
    """
    Compute correlations between interacting cell scores and module scores.

    Parameters
    ----------
    adata : AnnData
        Must contain:
        - ``uns['interacting_cell_results']`` (parametric or non-parametric interacting cell scores)
        - ``obsm['module_scores']`` or ``obsm['super_module_scores']``
        - ``uns['metabolites']`` or ``uns['gene_pairs_sig_names']``
    cor_method : {"pearson", "spearman"}, default "pearson"
        Statistical method used to compute correlations.
    test : {"parametric", "non-parametric"}
        Which interacting-cell score set to use.
        - `"parametric"` → uses ``uns['interacting_cell_results']['p']``
        - `"non-parametric"` → uses ``uns['interacting_cell_results']['np']``
    interaction_type : {"metabolite", "gene_pair"}, default "metabolite"
        Select whether to correlate:
        - metabolite scores, or
        - gene pair scores.
    only_sig_values : bool, default False
        If True, use only significant interacting cell score values (`cs_sig_pval` or `cs_sig_FDR`).
    normalize_values : bool, default False
        Apply min–max normalization to interacting cell score values per interaction.
    use_FDR : bool, default True
        If ``only_sig_values=True``, determines whether to filter by FDR or raw p-values.
    use_super_modules : bool, default False
        Whether to use super-module scores (``obsm['super_module_scores']``) instead of module scores.
    """
    MODULE_KEY = "super_module_scores" if use_super_modules else "module_scores"

    if cor_method not in ["pearson", "spearman"]:
        raise ValueError(f'Invalid method: {cor_method}. Choose either "pearson" or "spearman".')

    adata.uns["cor_method"] = cor_method

    if test not in ["parametric", "non-parametric"]:
        raise ValueError('The "test" variable should be one of ["parametric", "non-parametric"].')

    test_str = "p" if test == "parametric" else "np"

    if interaction_type not in ["metabolite", "gene_pair"]:
        raise ValueError(
            'The "interaction_type" variable should be one of ["metabolite", "gene_pair"].'
        )

    interaction_type_str = "m" if interaction_type == "metabolite" else "gp"

    if only_sig_values:
        sig_str = "FDR" if use_FDR else "pval"
        interaction_scores = adata.uns["interacting_cell_results"][test_str][interaction_type_str][
            f"cs_sig_{sig_str}"
        ]
    else:
        interaction_scores = adata.uns["interacting_cell_results"][test_str][interaction_type_str][
            "cs"
        ]

    if normalize_values:
        interaction_scores = interaction_scores.apply(
            lambda x: (x - x.min()) / (x.max() - x.min()), axis=0
        )  # We apply min-max normalization

    interaction_type_names_key = (
        "metabolites" if interaction_type == "metabolite" else "gene_pairs_sig_names"
    )
    interaction_scores = pd.DataFrame(
        interaction_scores, index=adata.obs_names, columns=adata.uns[interaction_type_names_key]
    )

    metabolites = interaction_scores.columns.tolist()
    modules = adata.obsm[MODULE_KEY].columns.tolist()

    cor_pval_df = pd.DataFrame(index=modules)
    cor_coef_df = pd.DataFrame(index=modules)

    for metab in metabolites:
        correlation_values = []
        pvals = []

        for module in modules:
            metab_df = interaction_scores[metab]
            module_df = adata.obsm[MODULE_KEY][module]

            if cor_method == "pearson":
                correlation_value, pval = pearsonr(metab_df, module_df)
            elif cor_method == "spearman":
                correlation_value, pval = spearmanr(metab_df, module_df)

            correlation_values.append(correlation_value)
            pvals.append(pval)

        cor_coef_df[metab] = correlation_values
        cor_pval_df[metab] = pvals

    cor_pval_df = cor_pval_df.replace(np.nan, 1)
    cor_coef_df = cor_coef_df.replace(np.nan, 0)
    cor_FDR_values = multipletests(cor_pval_df.values.flatten(), method="fdr_bh")[1]
    cor_FDR_df = pd.DataFrame(
        cor_FDR_values.reshape(cor_pval_df.shape),
        index=cor_pval_df.index,
        columns=cor_pval_df.columns,
    )

    adata.uns["interaction_module_correlation_coefs"] = cor_coef_df
    adata.uns["interaction_module_correlation_pvals"] = cor_pval_df
    adata.uns["interaction_module_correlation_FDR"] = cor_FDR_df

    return
