import numpy as np
import logging
import pandas as pd
import anndata
import scipy.sparse as sp_sparse
import torch

from typing import Optional

from ._anndata_utils import _check_nonnegative_integers

logger = logging.getLogger(__name__)


def highly_variable_genes_seurat_v3(
    adata: anndata.AnnData,
    layer: Optional[str] = None,
    n_top_genes: int = 2000,
    batch_key: Optional[str] = None,
    span: Optional[float] = 0.3,
    subset: bool = False,
    inplace: bool = True,
) -> Optional[pd.DataFrame]:
    """\
    Implements highly variable gene selection using the `vst` method in Seurat v3.

    Expects count data.

    This is a temporary implementation that will be replaced by the one in Scanpy.
    For further implemenation details see https://www.overleaf.com/read/ckptrbgzzzpg

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    layer
        If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.
    n_top_genes
        Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'`.
    span
        The fraction of the data (cells) used when estimating the variance in the loess
        model fit if `flavor='seurat_v3'`.
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    batch_key
        If specified, highly-variable genes are selected within each batch separately and merged.
        This simple process avoids the selection of batch-specific genes and acts as a
        lightweight batch correction method. Genes are first sorted by how many batches
        they are a HVG. Ties are broken by the median (across batches) rank based on
        within-batch normalized variance.

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pd.DataFrame`) or
    updates `.var` with the following fields

    highly_variable : bool
        boolean indicator of highly-variable genes
    **means**
        means per gene
    **variances**
        variance per gene
    **variances_norm**
        normalized variance per gene, averaged in the case of multiple batches
    highly_variable_rank : float
        Rank of the gene according to normalized variance, median rank in the case of multiple batches
    highly_variable_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as HVG
    """
    from scanpy.preprocessing._utils import _get_mean_var

    try:
        from skmisc.loess import loess
    except ImportError:
        raise ImportError(
            "Please install skmisc package via `pip install --user scikit-misc"
        )

    X = adata.layers[layer] if layer is not None else adata.X
    if _check_nonnegative_integers(X) is False:
        raise ValueError(
            "`pp.highly_variable_genes` with `flavor='seurat_v3'` expects "
            "raw count data."
        )

    if batch_key is None:
        batch_info = pd.Categorical(np.zeros(adata.shape[0], dtype=int))
    else:
        batch_info = adata.obs[batch_key].values

    norm_gene_vars = []
    for b in np.unique(batch_info):

        ad = adata[batch_info == b]
        X = ad.layers[layer] if layer is not None else ad.X

        mean, var = _get_mean_var(X)
        not_const = var > 0
        estimat_var = np.zeros(adata.shape[1], dtype=np.float64)

        y = np.log10(var[not_const])
        x = np.log10(mean[not_const])
        model = loess(x, y, span=span, degree=2)
        model.fit()
        estimat_var[not_const] = model.outputs.fitted_values
        reg_std = np.sqrt(10 ** estimat_var)

        batch_counts = X.astype(np.float64).copy()
        # clip large values as in Seurat
        N = np.sum(batch_info == b)
        vmax = np.sqrt(N)
        clip_val = reg_std * vmax + mean
        if sp_sparse.issparse(batch_counts):
            batch_counts = sp_sparse.csr_matrix(batch_counts)
            mask = batch_counts.data > clip_val[batch_counts.indices]
            batch_counts.data[mask] = clip_val[batch_counts.indices[mask]]
        else:
            clip_val_broad = np.broadcast_to(clip_val, batch_counts.shape)
            np.putmask(batch_counts, batch_counts > clip_val_broad, clip_val_broad)

        if sp_sparse.issparse(batch_counts):
            squared_batch_counts_sum = np.array(batch_counts.power(2).sum(axis=0))
            batch_counts_sum = np.array(batch_counts.sum(axis=0))
        else:
            squared_batch_counts_sum = np.square(batch_counts).sum(axis=0)
            batch_counts_sum = batch_counts.sum(axis=0)

        norm_gene_var = (1 / ((N - 1) * np.square(reg_std))) * (
            (N * np.square(mean))
            + squared_batch_counts_sum
            - 2 * batch_counts_sum * mean
        )
        norm_gene_vars.append(norm_gene_var.reshape(1, -1))

    norm_gene_vars = np.concatenate(norm_gene_vars, axis=0)
    # argsort twice gives ranks, small rank means most variable
    ranked_norm_gene_vars = np.argsort(np.argsort(-norm_gene_vars, axis=1), axis=1)

    # this is done in SelectIntegrationFeatures() in Seurat v3
    ranked_norm_gene_vars = ranked_norm_gene_vars.astype(np.float32)
    num_batches_high_var = np.sum(
        (ranked_norm_gene_vars < n_top_genes).astype(int), axis=0
    )
    ranked_norm_gene_vars[ranked_norm_gene_vars >= n_top_genes] = np.nan
    ma_ranked = np.ma.masked_invalid(ranked_norm_gene_vars)
    median_ranked = np.ma.median(ma_ranked, axis=0).filled(np.nan)

    df = pd.DataFrame(index=np.array(adata.var_names))
    df["highly_variable_nbatches"] = num_batches_high_var
    df["highly_variable_rank"] = median_ranked
    df["variances_norm"] = np.mean(norm_gene_vars, axis=0)
    df["means"] = mean
    df["variances"] = var

    df.sort_values(
        ["highly_variable_rank", "highly_variable_nbatches"],
        ascending=[True, False],
        na_position="last",
        inplace=True,
    )
    df["highly_variable"] = False
    df.loc[: int(n_top_genes), "highly_variable"] = True
    df = df.loc[adata.var_names]

    if inplace or subset:
        adata.uns["hvg"] = {"flavor": "seurat_v3"}
        logger.info(
            "added\n"
            "    'highly_variable', boolean vector (adata.var)\n"
            "    'highly_variable_rank', float vector (adata.var)\n"
            "    'means', float vector (adata.var)\n"
            "    'variances', float vector (adata.var)\n"
            "    'variances_norm', float vector (adata.var)"
        )
        adata.var["highly_variable"] = df["highly_variable"].values
        adata.var["highly_variable_rank"] = df["highly_variable_rank"].values
        adata.var["means"] = df["means"].values
        adata.var["variances"] = df["variances"].values
        adata.var["variances_norm"] = df["variances_norm"].values.astype(
            "float64", copy=False
        )
        if batch_key is not None:
            adata.var["highly_variable_nbatches"] = df[
                "highly_variable_nbatches"
            ].values
        if subset:
            adata._inplace_subset_var(df["highly_variable"].values)
    else:
        if batch_key is None:
            df = df.drop(["highly_variable_nbatches"], axis=1)
        return df


def poisson_gene_selection(
    adata,
    layer: Optional[str] = None,
    n_top_genes: int = 4000,
    use_cuda: bool = True,
    subset: bool = False,
    inplace: bool = True,
    n_samples: int = 10000,
    batch_key: str = None,
    silent: bool = False,
    minibatch_size: int = 5000,
    **kwargs,
):
    """\
    Rank and select genes based on the enrichment of zero counts in data compared to a Poisson count model.

    This is based on M3Drop: https://github.com/tallulandrews/M3Drop

    The method accounts for library size internally, a raw count matrix should be provided.

    Instead of Z-test, enrichment of zeros is quantified by posterior
    probabilites from a binomial model, computed through sampling.


    Parameters
    ----------
    adata
        AnnData object (with sparse X matrix).
    layer
        If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.
    n_top_genes
        How many variable genes to select.
    use_cuda
        Whether to use GPU
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    n_samples
        The number of Binomial samples to use to estimate posterior probability
        of enrichment of zeros for each gene.
    batch_key
        key in adata.obs that contains batch info. If None, do not use batch info.
        Defatult: ``None``.
    silent
        If ``True``, disables the progress bar.
    minibatch_size
        Size of temporary matrix for incremental calculation. Larger is faster but
        requires more RAM or GPU memory. (The default should be fine unless
        there are hundreds of millions cells or millions of genes.)

    Returns
    -------
    Depending on `inplace` returns calculated metrics (:class:`~pd.DataFrame`) or
    updates `.var` with the following fields

    highly_variable : bool
        boolean indicator of highly-variable genes
    **observed_fraction_zeros**
        fraction of observed zeros per gene
    **expected_fraction_zeros**
        expected fraction of observed zeros per gene
    prob_zero_enrichment : float
        Probability of zero enrichment, median across batches in the case of multiple batches
    prob_zero_enrichment_rank : float
        Rank of the gene according to probability of zero enrichment, median rank in the case of multiple batches
    prob_zero_enriched_nbatches : int
        If batch_key is given, this denotes in how many batches genes are detected as zero enriched

    """
    from tqdm import tqdm
    import sys

    X = adata.layers[layer] if layer is not None else adata.X
    if _check_nonnegative_integers(X) is False:
        raise ValueError("`poisson_gene_selection` expects " "raw count data.")

    use_cuda = use_cuda and torch.cuda.is_available()
    dev = torch.device("cuda:0" if use_cuda else "cpu")

    if batch_key is None:
        batch_info = pd.Categorical(np.zeros(adata.shape[0], dtype=int))
    else:
        batch_info = adata.obs[batch_key]

    prob_zero_enrichments = []
    obs_frac_zeross = []
    exp_frac_zeross = []
    for b in np.unique(batch_info):

        ad = adata[batch_info == b]
        X = ad.layers[layer] if layer is not None else ad.X

        # Calculate empirical statistics.
        scaled_means = torch.from_numpy(X.sum(0) / X.sum())[0].to(dev)
        total_counts = torch.from_numpy(X.sum(1))[:, 0].to(dev)

        observed_fraction_zeros = torch.from_numpy(1.0 - (X > 0).sum(0) / X.shape[0])[
            0
        ].to(dev)

        # Calculate probability of zero for a Poisson model.
        # Perform in batches to save memory.
        minibatch_size = min(total_counts.shape[0], minibatch_size)
        n_batches = total_counts.shape[0] // minibatch_size

        expected_fraction_zeros = torch.zeros(scaled_means.shape).to(dev)

        for i in range(n_batches):
            total_counts_batch = total_counts[
                i * minibatch_size : (i + 1) * minibatch_size
            ]
            # Use einsum for outer product.
            expected_fraction_zeros += torch.exp(
                -torch.einsum("i,j->ij", [scaled_means, total_counts_batch])
            ).sum(1)

        total_counts_batch = total_counts[(i + 1) * minibatch_size :]
        expected_fraction_zeros += torch.exp(
            -torch.einsum("i,j->ij", [scaled_means, total_counts_batch])
        ).sum(1)
        expected_fraction_zeros /= X.shape[0]

        # Compute probability of enriched zeros through sampling from Binomial distributions.
        observed_zero = torch.distributions.Binomial(probs=observed_fraction_zeros)
        expected_zero = torch.distributions.Binomial(probs=expected_fraction_zeros)

        extra_zeros = torch.zeros(expected_fraction_zeros.shape).to(dev)
        for i in tqdm(
            range(n_samples),
            disable=silent,
            desc="Sampling from binomial",
            file=sys.stdout,
        ):
            extra_zeros += observed_zero.sample() > expected_zero.sample()

        prob_zero_enrichment = (extra_zeros / n_samples).cpu().numpy()

        obs_frac_zeros = observed_fraction_zeros.cpu().numpy()
        exp_frac_zeros = expected_fraction_zeros.cpu().numpy()

        # Clean up memory (tensors seem to stay in GPU unless actively deleted).
        del scaled_means
        del total_counts
        del expected_fraction_zeros
        del observed_fraction_zeros
        del extra_zeros

        if use_cuda:
            torch.cuda.empty_cache()

        prob_zero_enrichments.append(prob_zero_enrichment.reshape(1, -1))
        obs_frac_zeross.append(obs_frac_zeros.reshape(1, -1))
        exp_frac_zeross.append(exp_frac_zeros.reshape(1, -1))

    # Combine per batch results

    prob_zero_enrichments = np.concatenate(prob_zero_enrichments, axis=0)
    obs_frac_zeross = np.concatenate(obs_frac_zeross, axis=0)
    exp_frac_zeross = np.concatenate(exp_frac_zeross, axis=0)

    ranked_prob_zero_enrichments = prob_zero_enrichments.argsort(axis=1).argsort(axis=1)
    median_prob_zero_enrichments = np.median(prob_zero_enrichments, axis=0)

    median_obs_frac_zeross = np.median(obs_frac_zeross, axis=0)
    median_exp_frac_zeross = np.median(exp_frac_zeross, axis=0)

    median_ranked = np.median(ranked_prob_zero_enrichments, axis=0)

    num_batches_zero_enriched = np.sum(
        ranked_prob_zero_enrichments >= (adata.shape[1] - n_top_genes), axis=0
    )

    df = pd.DataFrame(index=np.array(adata.var_names))
    df["observed_fraction_zeros"] = median_obs_frac_zeross
    df["expected_fraction_zeros"] = median_exp_frac_zeross
    df["prob_zero_enriched_nbatches"] = num_batches_zero_enriched
    df["prob_zero_enrichment"] = median_prob_zero_enrichments
    df["prob_zero_enrichment_rank"] = median_ranked

    df["highly_variable"] = False
    sort_columns = ["prob_zero_enriched_nbatches", "prob_zero_enrichment_rank"]
    top_genes = df.nlargest(n_top_genes, sort_columns).index
    df.loc[top_genes, "highly_variable"] = True

    if inplace or subset:
        adata.uns["hvg"] = {"flavor": "poisson_zeros"}
        logger.debug(
            "added\n"
            "    'highly_variable', boolean vector (adata.var)\n"
            "    'prob_zero_enrichment_rank', float vector (adata.var)\n"
            "    'prob_zero_enrichment' float vector (adata.var)\n"
            "    'observed_fraction_zeros', float vector (adata.var)\n"
            "    'expected_fraction_zeros', float vector (adata.var)\n"
        )
        adata.var["highly_variable"] = df["highly_variable"].values
        adata.var["observed_fraction_zeros"] = df["observed_fraction_zeros"].values
        adata.var["expected_fraction_zeros"] = df["expected_fraction_zeros"].values
        adata.var["prob_zero_enriched_nbatches"] = df[
            "prob_zero_enriched_nbatches"
        ].values
        adata.var["prob_zero_enrichment"] = df["prob_zero_enrichment"].values
        adata.var["prob_zero_enrichment_rank"] = df["prob_zero_enrichment_rank"].values

        if batch_key is not None:
            adata.var["prob_zero_enriched_nbatches"] = df[
                "prob_zero_enriched_nbatches"
            ].values
        if subset:
            adata._inplace_subset_var(df["highly_variable"].values)
    else:
        if batch_key is None:
            df = df.drop(["prob_zero_enriched_nbatches"], axis=1)
        return df


def organize_cite_seq_10x(adata: anndata.AnnData, copy: bool = False):
    """\
    Organize anndata object loaded from 10x for scvi models

    Parameters
    ----------
    adata
        AnnData object with RNA and protein data in `.X`
    copy
        Whether to copy the anndata object

    Returns
    -------
    If copy is True, returns anndata object organized for input to scvi models

    Else, updates the anndata inplace

    """

    if copy:
        adata = adata.copy()

    pro_array = adata[:, adata.var["feature_types"] == "Antibody Capture"].X.A
    pro_names = np.array(
        adata.var_names[adata.var["feature_types"] == "Antibody Capture"]
    )
    pro_df = pd.DataFrame(pro_array, index=adata.obs_names, columns=pro_names)
    adata.obsm["protein_expression"] = pro_df

    genes = (adata.var["feature_types"] != "Antibody Capture").values
    adata._inplace_subset_var(genes)

    if copy:
        return adata
