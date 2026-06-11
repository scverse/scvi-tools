import time

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm, zscore
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from scvi.external.harreman._utils import _resolve_device
from scvi.external.harreman.preprocessing.anndata import counts_from_anndata
from scvi.external.harreman.tools.knn import make_weights_non_redundant

from .local_autocorrelation import center_counts_torch


def calculate_module_scores(
    adata: AnnData,
    device: torch.device | str = "auto",
    verbose: bool | None = False,
) -> None:
    """
    Calculate module scores for gene modules across cells.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData). Required fields in `adata.uns`:
        - 'layer_key': name of the layer from which to extract the expression matrix
        - 'model': statistical model used to normalize expression (e.g., 'DANB', 'normal')
        - 'umi_counts': total UMI counts per cell
        - 'gene_modules_dict': dictionary mapping module IDs (as strings) to lists of gene names
    device : torch.device, optional
        Device to use for computation (e.g., CUDA or CPU). Defaults to GPU if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        The following results are stored in the `AnnData` object:
        - `adata.obsm['module_scores']`: (cells x modules) DataFrame with per-cell module scores
        - `adata.varm['gene_loadings']`: (genes x modules) DataFrame with gene loadings per module
        - `adata.uns['gene_modules']`: dictionary mapping module names to gene lists
    """
    start = time.time()
    device = _resolve_device(device)
    device = _resolve_device(device)

    layer_key = adata.uns["layer_key"]
    model = adata.uns["model"]

    use_raw = layer_key == "use_raw"
    modules = adata.uns["gene_modules_dict"].copy()

    umi_counts = adata.uns["umi_counts"]

    modules_to_compute = sorted([x for x in modules.keys() if x != "-1"])
    mod_list = [int(mod) for mod in modules_to_compute]
    mod_list.sort()
    modules_to_compute = [str(mod) for mod in mod_list]

    if verbose:
        print(f"Computing scores for {len(modules_to_compute)} modules...")

    module_scores = {}
    gene_loadings = pd.DataFrame(index=adata.var_names)
    gene_modules = {}
    for module in tqdm(modules_to_compute):
        module_genes = modules[module]

        scores, loadings = compute_scores(
            adata[:, module_genes],
            layer_key,
            model,
            umi_counts,
            device,
        )

        module_name = f"Module {module}" if "Module" not in module else module
        module_scores[module_name] = scores
        gene_loadings[module_name] = pd.Series(loadings, index=module_genes)
        gene_modules[module_name] = module_genes

    module_scores = pd.DataFrame(module_scores)

    module_scores.index = adata.obs_names if not use_raw else adata.raw.obs.index

    adata.obsm["module_scores"] = module_scores
    adata.varm["gene_loadings"] = gene_loadings
    adata.uns["gene_modules"] = gene_modules

    if verbose:
        print("Finished computing module scores in %.3f seconds" % (time.time() - start))

    return


def compute_scores(adata, layer_key, model, num_umi, device, _lambda=0.9):
    """counts_sub: row-subset of counts matrix with genes in the module"""
    # Get the weights matrix
    weights = make_weights_non_redundant(adata.obsp["weights"]).tocoo()
    weights = torch.sparse_coo_tensor(
        torch.tensor(np.vstack((weights.row, weights.col)), dtype=torch.long, device=device),
        torch.tensor(weights.data, dtype=torch.float64, device=device),
        torch.Size(weights.shape),
        device=device,
    )

    # Get gene expression counts for the module (dense)
    counts_sub = counts_from_anndata(adata, layer_key, dense=True)

    # Convert to tensors
    num_umi = torch.tensor(num_umi, dtype=torch.float64, device=device)
    counts_sub = torch.tensor(counts_sub, dtype=torch.float64, device=device)

    # Center values
    sample_specific = "sample_key" in adata.uns.keys()
    if sample_specific:
        sample_key = adata.uns["sample_key"]
        for sample in adata.obs[sample_key].unique():
            subset = np.where(adata.obs[sample_key] == sample)[0]
            counts_sub[:, subset] = center_counts_torch(
                counts_sub[:, subset], num_umi[subset], model
            )
    else:
        counts_sub = center_counts_torch(counts_sub, num_umi, model)

    # Smooth the counts using weights
    out = torch.matmul(
        weights + weights.transpose(0, 1), counts_sub.T
    )  # (cells x cells) @ (cells x genes)^T = (cells x genes)^T
    weights_sum = (
        torch.sparse.sum(weights, dim=0).to_dense() + torch.sparse.sum(weights, dim=1).to_dense()
    )  # shape (cells,)
    weights_sum[weights_sum == 0] = 1.0
    out = out / weights_sum[:, None]  # normalize
    cc_smooth = _lambda * out.T + (1 - _lambda) * counts_sub  # (genes x cells)

    # Perform PCA on cells (transpose to cells x genes)
    pca_data = cc_smooth.T.cpu().numpy()
    pca = PCA(n_components=1)
    scores = pca.fit_transform(pca_data)
    loadings = pca.components_.T

    # Flip sign if needed
    if pca.components_.mean() < 0:
        scores *= -1
        loadings *= -1
    scores = scores[:, 0]
    loadings = loadings[:, 0]

    return scores, loadings


def sort_linkage(Z, node_index, node_values):
    """Sorts linkage by 'node_values' in place"""
    N = Z.shape[0] + 1  # number of leaves

    if node_index < 0:
        return

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    swap = False

    if left_child < 0 and right_child < 0:
        swap = False
    elif left_child < 0 and right_child >= 0:
        swap = True
    elif left_child >= 0 and right_child < 0:
        swap = False
    else:
        if node_values[left_child] > node_values[right_child]:
            swap = True
        else:
            swap = False

    if swap:
        Z[node_index, 0] = right_child + N
        Z[node_index, 1] = left_child + N

    sort_linkage(Z, left_child, node_values)
    sort_linkage(Z, right_child, node_values)


def calc_mean_dists(Z, node_index, out_mean_dists):
    """Calculate mean density of joins for sub-trees underneath each node."""
    N = Z.shape[0] + 1  # number of leaves

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        left_average = 0
        left_merges = 0
    else:
        left_average, left_merges = calc_mean_dists(Z, left_child, out_mean_dists)

    if right_child < 0:
        right_average = 0
        right_merges = 0
    else:
        right_average, right_merges = calc_mean_dists(Z, right_child, out_mean_dists)

    this_height = Z[node_index, 2]
    this_merges = left_merges + right_merges + 1
    this_average = (
        left_average * left_merges + right_average * right_merges + this_height
    ) / this_merges

    out_mean_dists[node_index] = this_average

    return this_average, this_merges


def prop_label(Z, node_index, label, labels, out_clusters):
    """Propagate node labels downward if they are not -1.

    Used to find the correct cluster label at the leaves.
    """
    N = Z.shape[0] + 1  # number of leaves

    if label == -1:
        label = labels[node_index]

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        out_clusters[left_child + N] = label
    else:
        prop_label(Z, left_child, label, labels, out_clusters)

    if right_child < 0:
        out_clusters[right_child + N] = label
    else:
        prop_label(Z, right_child, label, labels, out_clusters)


def prop_label2(Z, node_index, label, labels, out_clusters):
    """Propagate node labels downward.

    Helper method used in assign_modules.
    """
    N = Z.shape[0] + 1  # number of leaves

    parent_label = label
    this_label = labels[node_index]

    if this_label == -1:
        new_label = parent_label
    else:
        new_label = this_label

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        out_clusters[left_child + N] = new_label
    else:
        prop_label2(Z, left_child, new_label, labels, out_clusters)

    if right_child < 0:
        out_clusters[right_child + N] = new_label
    else:
        prop_label2(Z, right_child, new_label, labels, out_clusters)


def assign_modules(Z, leaf_labels, offset, MIN_THRESHOLD=10, Z_THRESHOLD=3):
    """Assign gene modules from a linkage matrix using mean-distance-guided cluster merging."""
    clust_i = 0

    labels = np.ones(Z.shape[0]) * -1
    N = Z.shape[0] + 1

    mean_dists = np.zeros(Z.shape[0])
    calc_mean_dists(Z, Z.shape[0] - 1, mean_dists)

    for i in range(Z.shape[0]):
        ca = int(Z[i, 0])
        cb = int(Z[i, 1])

        if ca - N < 0:  # leaf node
            n_members_a = 1
            clust_a = -1
        else:
            n_members_a = Z[ca - N, 3]
            clust_a = labels[ca - N]

        if cb - N < 0:  # leaf node
            n_members_b = 1
            clust_b = -1
        else:
            n_members_b = Z[cb - N, 3]
            clust_b = labels[cb - N]

        if Z[i, 2] > offset - Z_THRESHOLD:
            new_clust_assign = -1
        elif n_members_a >= MIN_THRESHOLD and n_members_b >= MIN_THRESHOLD:
            # don't join them
            # assign the one with the larger mean distance
            dist_a = mean_dists[ca - N]
            dist_b = mean_dists[cb - N]
            if dist_a >= dist_b:
                new_clust_assign = clust_a
            else:
                new_clust_assign = clust_b
        elif n_members_a >= MIN_THRESHOLD:
            new_clust_assign = clust_a
        elif n_members_b >= MIN_THRESHOLD:
            new_clust_assign = clust_b
        elif (n_members_b + n_members_a) >= MIN_THRESHOLD:
            # A new cluster is born!
            new_clust_assign = clust_i
            clust_i += 1
        else:
            new_clust_assign = -1  # Still too small

        labels[i] = new_clust_assign

    out_clusters = np.ones(N) * -2
    prop_label2(Z, Z.shape[0] - 1, labels[-1], labels, out_clusters)

    # remap out_clusters
    unique_clusters = list(np.sort(np.unique(out_clusters)))

    if -1 in unique_clusters:
        unique_clusters.remove(-1)

    clust_map = {x: i + 1 for i, x in enumerate(unique_clusters)}
    clust_map[-1] = -1

    out_clusters = [clust_map[x] for x in out_clusters]
    out_clusters = pd.Series(out_clusters, index=leaf_labels)

    return out_clusters


def assign_modules_core(Z, leaf_labels, offset, MIN_THRESHOLD=10, Z_THRESHOLD=3):
    """Assign core gene modules from a linkage matrix without mean-distance guidance."""
    clust_i = 0

    labels = np.ones(Z.shape[0]) * -1
    N = Z.shape[0] + 1

    for i in range(Z.shape[0]):
        ca = int(Z[i, 0])
        cb = int(Z[i, 1])

        if ca - N < 0:  # leaf node
            n_members_a = 1
            clust_a = -1
        else:
            n_members_a = Z[ca - N, 3]
            clust_a = labels[ca - N]

        if cb - N < 0:  # leaf node
            n_members_b = 1
            clust_b = -1
        else:
            n_members_b = Z[cb - N, 3]
            clust_b = labels[cb - N]

        if n_members_a >= MIN_THRESHOLD and n_members_b >= MIN_THRESHOLD:
            # don't join them
            new_clust_assign = -1
        elif Z[i, 2] > offset - Z_THRESHOLD:
            new_clust_assign = -1
        elif n_members_a >= MIN_THRESHOLD:
            new_clust_assign = clust_a
        elif n_members_b >= MIN_THRESHOLD:
            new_clust_assign = clust_b
        elif (n_members_b + n_members_a) >= MIN_THRESHOLD:
            # A new cluster is born!
            new_clust_assign = clust_i
            clust_i += 1
        else:
            new_clust_assign = -1  # Still too small

        labels[i] = new_clust_assign

    out_clusters = np.ones(N) * -2
    prop_label(Z, Z.shape[0] - 1, labels[-1], labels, out_clusters)

    # remap out_clusters
    unique_clusters = list(np.sort(np.unique(out_clusters)))

    if -1 in unique_clusters:
        unique_clusters.remove(-1)

    clust_map = {x: i + 1 for i, x in enumerate(unique_clusters)}
    clust_map[-1] = -1

    out_clusters = [clust_map[x] for x in out_clusters]
    out_clusters = pd.Series(out_clusters, index=leaf_labels)

    return out_clusters


def create_modules(
    adata: AnnData,
    min_gene_threshold: int | None = 20,
    fdr_threshold: float | None = 0.05,
    z_threshold: float | None = None,
    core_only: bool = False,
    verbose: bool | None = False,
):
    """
    Perform hierarchical clustering on gene-gene local correlation Z-scores to assign gene modules.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData) containing the local correlation Z-scores in
        `adata.uns['lc_zs']`.
    min_gene_threshold : int, optional (default: 20)
        Minimum number of genes required to define a module.
    fdr_threshold : float, optional (default: 0.05)
        FDR threshold used to determine the minimum Z-score significance if `z_threshold` is
        not provided.
    z_threshold : float, optional
        If provided, uses this Z-score as the cutoff for module inclusion instead of computing
        it from FDR.
    core_only : bool, optional (default: False)
        If True, assigns only tightly correlated (core) genes to modules and leaves others
        unassigned.

    Returns
    -------
    None
        The function modifies the `AnnData` object in place by adding the following to
        `adata.uns`:
        - `modules`: pandas Series mapping each gene to a module ID (integer, as string)
        - `gene_modules_dict`: dictionary mapping module IDs (as strings) to lists of gene names
        - `linkage`: linkage matrix from hierarchical clustering (for visualization or tree ops)
    """
    start = time.time()
    if verbose:
        print("Creating modules...")

    # Determine Z_Threshold from FDR threshold

    Z_scores = adata.uns["lc_zs"]

    if z_threshold is None:
        allZ = squareform(  # just in case slightly not symmetric
            Z_scores.values / 2 + Z_scores.values.T / 2
        )
        allZ = np.sort(allZ)
        allP = norm.sf(allZ)
        allP_c = multipletests(allP, method="fdr_bh")[1]
        ii = np.nonzero(allP_c < fdr_threshold)[0]
        if ii.size > 0:
            z_threshold = allZ[ii[0]]
        else:
            z_threshold = allZ[-1] + 1

    # Compute the linkage matrix
    dd = Z_scores.copy().values
    np.fill_diagonal(dd, 0)
    condensed = squareform(dd) * -1
    offset = condensed.min() * -1
    condensed += offset
    Z = linkage(condensed, method="average")

    # Linkage -> Modules
    if core_only:
        out_clusters = assign_modules_core(
            Z,
            offset=offset,
            MIN_THRESHOLD=min_gene_threshold,
            leaf_labels=Z_scores.index,
            Z_THRESHOLD=z_threshold,
        )
    else:
        out_clusters = assign_modules(
            Z,
            offset=offset,
            MIN_THRESHOLD=min_gene_threshold,
            leaf_labels=Z_scores.index,
            Z_THRESHOLD=z_threshold,
        )

    # Sort the leaves of the linkage matrix (for plotting)
    mean_dists = np.zeros(Z.shape[0])
    calc_mean_dists(Z, Z.shape[0] - 1, mean_dists)
    linkage_out = Z.copy()
    sort_linkage(linkage_out, Z.shape[0] - 1, mean_dists)

    out_clusters.name = "Module"

    gene_modules_dict = {}
    for mod in out_clusters.unique():
        gene_modules_dict[str(mod)] = out_clusters[out_clusters == mod].index.tolist()

    adata.uns["modules"] = out_clusters
    adata.uns["gene_modules_dict"] = gene_modules_dict
    adata.uns["linkage"] = linkage_out

    if verbose:
        print("Finished creating modules in %.3f seconds" % (time.time() - start))

    return


def compute_top_scoring_modules(
    adata: AnnData,
    sd: float | None = 1,
    use_super_modules: bool | None = False,
):
    """
    Identify the top-scoring module (or super-module) for each cell.

    Parameters
    ----------
    adata : AnnData
        Must contain a matrix of module or super-module scores in:
        - ``obsm['module_scores']`` or
        - ``obsm['super_module_scores']``
    sd : float, default 1
        Standard deviation threshold to determine strong module activation.
        Only modules with ``zscore > sd`` are considered strictly activated.
    use_super_modules : bool, default False
        If True, select among super-modules instead of standard modules.

    Returns
    -------
    pandas.Series
        A Series indexed by cell, containing the name of the top-scoring
        module/super-module for each cell.
    """
    MODULE_KEY = "super_module_scores" if use_super_modules else "module_scores"

    df = pd.DataFrame(
        zscore(adata.obsm[MODULE_KEY], axis=0),
        index=adata.obsm[MODULE_KEY].index,
        columns=adata.obsm[MODULE_KEY].columns,
    )

    top_scoring_modules = pd.Series(index=df.index)
    for mod_id, row in df.iterrows():
        above_threshold_low = row > 0
        above_threshold = row > sd
        if above_threshold.sum() == 1:
            top_scoring_modules[mod_id] = above_threshold.idxmax()
        else:
            highest_module = (
                row[above_threshold].idxmax()
                if above_threshold.sum() > 1
                else row.idxmax()
                if above_threshold_low.sum() > 0
                else np.nan
            )
            top_scoring_modules[mod_id] = highest_module

    return top_scoring_modules


def calculate_super_module_scores(
    adata: AnnData,
    super_module_dict: dict = None,
    device: torch.device | str = "auto",
    verbose: bool | None = False,
) -> None:
    """
    Calculate super-module scores for gene super-modules across cells.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (AnnData). Required fields in `adata.uns`:
        - 'layer_key': name of the layer from which to extract the expression matrix
        - 'model': statistical model used to normalize expression (e.g., 'DANB', 'normal')
        - 'umi_counts': total UMI counts per cell
        - 'gene_modules': dictionary mapping module IDs (as strings) to lists of gene names
    super_module_dict: dict
        Dictionary with super-module IDs (integers) as keys and lists of associated modules
        (as integers) as values.
    device : torch.device, optional
        Device to use for computation (e.g., CUDA or CPU). Defaults to GPU if available.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    Returns
    -------
    None
        The following results are stored in the `AnnData` object:
        - `adata.obsm['super_module_scores']`: (cells x super-modules) DataFrame with
          per-cell super-module activity scores
        - `adata.varm['gene_loadings_sm']`: (genes x super-modules) DataFrame with gene
          loadings for each super-module
        - `adata.uns['gene_modules_sm']`: dictionary mapping super-module names to gene lists
    """
    start = time.time()
    device = _resolve_device(device)

    gene_modules = adata.uns["gene_modules"]

    reverse_mapping = {value: key for key, values in super_module_dict.items() for value in values}
    adata.uns["super_modules"] = adata.uns["modules"].replace(reverse_mapping)

    super_module_dict = {key: values for key, values in super_module_dict.items() if key != -1}

    layer_key = adata.uns["layer_key"]
    model = adata.uns["model"]

    use_raw = layer_key == "use_raw"

    umi_counts = adata.uns["umi_counts"]

    if verbose:
        print(f"Computing scores for {len(super_module_dict.keys())} super-modules...")

    super_module_scores = {}
    gene_loadings_sm = pd.DataFrame(index=adata.var_names)
    gene_modules_sm = {}
    for sm, modules in tqdm(super_module_dict.items()):
        super_module = f"Module {sm}"
        modules = [f"Module {str(mod)}" for mod in modules]
        super_module_genes = [item for key in modules for item in gene_modules.get(key, [])]

        if len(super_module_genes) == 0:
            missing = [k for k in modules if k not in gene_modules]
            print(
                f"Warning: {super_module} skipped — no genes found. "
                f"Modules not in gene_modules: {missing}. "
                f"Available modules: {sorted(gene_modules.keys())}"
            )
            continue

        scores, loadings = compute_scores(
            adata[:, super_module_genes],
            layer_key,
            model,
            umi_counts,
            device,
        )

        super_module_scores[super_module] = scores
        gene_loadings_sm[super_module] = pd.Series(loadings, index=super_module_genes)
        gene_modules_sm[super_module] = super_module_genes

    super_module_scores = pd.DataFrame(super_module_scores)
    if not super_module_scores.empty:
        super_module_scores.index = adata.obs_names if not use_raw else adata.raw.obs.index

    adata.varm["gene_loadings_sm"] = gene_loadings_sm
    adata.obsm["super_module_scores"] = super_module_scores
    adata.uns["gene_modules_sm"] = gene_modules_sm

    print("Finished computing super-module scores in %.3f seconds" % (time.time() - start))

    return
