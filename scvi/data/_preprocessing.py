import logging
import tempfile
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import torch

from scvi._decorators import dependencies
from scvi.utils import track

from ._utils import _check_nonnegative_integers

logger = logging.getLogger(__name__)


@torch.inference_mode()
def poisson_gene_selection(
    adata,
    layer: Optional[str] = None,
    n_top_genes: int = 4000,
    use_gpu: bool = True,
    subset: bool = False,
    inplace: bool = True,
    n_samples: int = 10000,
    batch_key: str = None,
    silent: bool = False,
    minibatch_size: int = 5000,
    **kwargs,
) -> Optional[pd.DataFrame]:
    """
    Rank and select genes based on the enrichment of zero counts.

    Enrichment is considered by comparing data to a Poisson count model.
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
    use_gpu
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
    data = adata.layers[layer] if layer is not None else adata.X
    if _check_nonnegative_integers(data) is False:
        raise ValueError("`poisson_gene_selection` expects " "raw count data.")

    use_gpu = use_gpu and torch.cuda.is_available()

    if batch_key is None:
        batch_info = pd.Categorical(np.zeros(adata.shape[0], dtype=int))
    else:
        batch_info = adata.obs[batch_key]

    prob_zero_enrichments = []
    obs_frac_zeross = []
    exp_frac_zeross = []
    for b in np.unique(batch_info):

        ad = adata[batch_info == b]
        data = ad.layers[layer] if layer is not None else ad.X

        # Calculate empirical statistics.
        sum_0 = np.asarray(data.sum(0)).ravel()
        scaled_means = torch.from_numpy(sum_0 / sum_0.sum())
        if use_gpu is True:
            scaled_means = scaled_means.cuda()
        dev = scaled_means.device
        total_counts = torch.from_numpy(np.asarray(data.sum(1)).ravel()).to(dev)

        observed_fraction_zeros = torch.from_numpy(
            np.asarray(1.0 - (data > 0).sum(0) / data.shape[0]).ravel()
        ).to(dev)

        # Calculate probability of zero for a Poisson model.
        # Perform in batches to save memory.
        minibatch_size = min(total_counts.shape[0], minibatch_size)
        n_batches = total_counts.shape[0] // minibatch_size

        expected_fraction_zeros = torch.zeros(scaled_means.shape, device=dev)

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
        expected_fraction_zeros /= data.shape[0]

        # Compute probability of enriched zeros through sampling from Binomial distributions.
        observed_zero = torch.distributions.Binomial(probs=observed_fraction_zeros)
        expected_zero = torch.distributions.Binomial(probs=expected_fraction_zeros)

        extra_zeros = torch.zeros(expected_fraction_zeros.shape, device=dev)
        for _ in track(
            range(n_samples),
            description="Sampling from binomial...",
            disable=silent,
            style="tqdm",  # do not change
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

        if use_gpu:
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


def organize_cite_seq_10x(
    adata: anndata.AnnData, copy: bool = False
) -> Optional[anndata.AnnData]:
    """
    Organize anndata object loaded from 10x for scvi models.

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

    Examples
    --------
    >>> adata = scanpy.read_10x_h5(<path_to_10x_h5_file>, gex_only=False)
    >>> adata
    AnnData object with n_obs × n_vars = 713 × 33555
        var: 'gene_ids', 'feature_types', 'genome'
    >>> organize_cite_seq_10x(adata)
    >>> adata
    AnnData object with n_obs × n_vars = 713 × 33538
        var: 'gene_ids', 'feature_types', 'genome'
        obsm: 'protein_expression'
    """
    if copy:
        adata = adata.copy()

    pro_array = adata[:, adata.var["feature_types"] == "Antibody Capture"].X.copy().A
    pro_names = np.array(
        adata.var_names[adata.var["feature_types"] == "Antibody Capture"]
    )

    genes = (adata.var["feature_types"] != "Antibody Capture").values
    adata._inplace_subset_var(genes)

    pro_df = pd.DataFrame(pro_array, index=adata.obs_names, columns=pro_names)
    adata.obsm["protein_expression"] = pro_df

    if copy:
        return adata


def organize_multiome_anndatas(
    multi_anndata: anndata.AnnData,
    rna_anndata: Optional[anndata.AnnData] = None,
    atac_anndata: Optional[anndata.AnnData] = None,
    modality_key: str = "modality",
) -> anndata.AnnData:
    """
    Concatenate multiome and single-modality input anndata objects.

    These anndata objects should already have been preprocessed so that both single-modality
    objects use a subset of the features used in the multiome object. The feature names (index of
    `.var`) should match between the objects.

    Parameters
    ----------
    multi_anndata
        AnnData object with Multiome data (Gene Expression and Chromatin Accessibility)
    rna_anndata
        AnnData object with gene expression data
    atac_anndata
        AnnData object with chromatin accessibility data
    modality_key
        The key to add to the resulting AnnData `.obs`, indicating the modality each cell originated
        from. Default is "modality".

    Notes
    -----
    Features that exist in either rna_anndata or atac_anndata but do not exist in multi_anndata will
    be discarded.

    Returns
    -------
    An AnnData object with all cells in the input objects
    """
    res_anndata = multi_anndata.copy()

    modality_ann = ["paired"] * multi_anndata.shape[0]
    obs_names = list(multi_anndata.obs.index.values)

    def _concat_anndata(multi_anndata, other):
        shared_features = np.intersect1d(
            other.var.index.values, multi_anndata.var.index.values
        )
        if not len(shared_features) > 0:
            raise ValueError("No shared features between Multiome and other AnnData.")

        other = other[:, shared_features]
        return multi_anndata.concatenate(other, join="outer", batch_key=modality_key)

    if rna_anndata is not None:
        res_anndata = _concat_anndata(res_anndata, rna_anndata)

        modality_ann += ["expression"] * rna_anndata.shape[0]
        obs_names += list(rna_anndata.obs.index.values)

    if atac_anndata is not None:
        res_anndata = _concat_anndata(res_anndata, atac_anndata)

        modality_ann += ["accessibility"] * atac_anndata.shape[0]
        obs_names += list(atac_anndata.obs.index.values)

    # set .obs stuff
    res_anndata.obs[modality_key] = modality_ann
    res_anndata.obs.index = (
        pd.Series(obs_names) + "_" + res_anndata.obs[modality_key].values
    )

    # keep the feature order as the original order in the multiomic anndata
    res_anndata = res_anndata[:, multi_anndata.var.index.values]
    res_anndata.var = multi_anndata.var.copy()
    return res_anndata.copy()


def _dna_to_code(nt: str) -> int:
    nt = nt.upper()
    if nt == "A":
        return 0
    elif nt == "C":
        return 1
    elif nt == "G":
        return 2
    elif nt == "T":
        return 3
    else:
        # scBasset does this
        return np.random.randint(0, 3)


@dependencies("genomepy")
def add_dna_sequence(
    adata: anndata.AnnData,
    seq_len: int = 1334,
    genome_name: str = "hg38",
    genome_dir: Optional[Path] = None,
    genome_provider: str = "UCSC",
    install_genome: bool = True,
    chr_var_key: str = "chr",
    start_var_key: str = "start",
    end_var_key: str = "end",
    sequence_varm_key: str = "dna_sequence",
    code_varm_key: str = "dna_code",
) -> None:
    """
    Add DNA sequence to AnnData object.

    Uses genomepy under the hood to download the genome.

    Parameters
    ----------
    adata
        AnnData object with chromatin accessiblity data
    seq_len
        Length of DNA sequence to extract around peak center.
        Defaults to value used in scBasset.
    genome_name
        Name of genome to use, installed with genomepy
    genome_dir
        Directory to install genome to, if not already installed
    genome_provider
        Provider of genome, passed to genomepy
    install_genome
        Install the genome with genomepy. If False, `genome_provider` is not used,
        and a genome is loaded with `genomepy.Genome(genome_name, genomes_dir=genome_dir)`
    chr_var_key
        Key in `.var` for chromosome
    start_var_key
        Key in `.var` for start position
    end_var_key
        Key in `.var` for end position
    sequence_varm_key
        Key in `.varm` for added DNA sequence
    code_varm_key
        Key in `.varm` for added DNA sequence, encoded as integers

    Returns
    -------
    None

    Adds fields to `.varm`:
        sequence_varm_key: DNA sequence
        code_varm_key: DNA sequence, encoded as integers
    """
    import genomepy

    if genome_dir is None:
        tempdir = tempfile.TemporaryDirectory()
        genome_dir = tempdir.name

    if install_genome:
        g = genomepy.install_genome(
            genome_name, genome_provider, genomes_dir=genome_dir
        )
    else:
        g = genomepy.Genome(genome_name, genomes_dir=genome_dir)

    output_dfs = []
    chroms = adata.var[chr_var_key].unique()
    df = adata.var[[chr_var_key, start_var_key, end_var_key]]
    for chrom in track(chroms):
        chrom_df = df[df[chr_var_key] == chrom]
        lengths = chrom_df[end_var_key] - chrom_df[start_var_key]

        block_starts = (
            chrom_df[start_var_key] + (lengths // 2.0) - (seq_len // 2.0)
        ).astype(int)
        block_ends = block_starts + seq_len

        concat_seq = str(
            g.get_spliced_seq(chrom, zip(block_starts, block_ends - 1))
        ).upper()
        concat_seq = [iter(concat_seq)] * seq_len
        concat_seq = list(zip(*concat_seq))
        assert len(concat_seq) == len(chrom_df)

        output_dfs.append(pd.DataFrame(np.array(concat_seq), index=chrom_df.index))

    output_df = pd.concat(output_dfs, axis=0).loc[adata.var_names]
    adata.varm[sequence_varm_key] = output_df
    adata.varm[code_varm_key] = output_df.applymap(_dna_to_code)
