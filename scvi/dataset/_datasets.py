import anndata

from typing import Optional, List
from scvi.dataset._built_in_data._brain_large import _load_brainlarge_dataset
from scvi.dataset._built_in_data._cortex import _load_cortex
from scvi.dataset._built_in_data._csv import (
    _load_breast_cancer_dataset,
    _load_mouse_ob_dataset,
)
from scvi.dataset._built_in_data._synthetic import _generate_synthetic
from scvi.dataset._built_in_data._cite_seq import _load_pbmcs_10x_cite_seq
from scvi.dataset._built_in_data._loom import (
    _load_annotation_simulation,
    _load_frontalcortex_dropseq,
    _load_prefrontalcortex_starmap,
    _load_retina,
)
from scvi.dataset._built_in_data._pbmc import (
    _load_purified_pbmc_dataset,
    _load_pbmc_dataset,
)
from scvi.dataset._built_in_data._seqfish import _load_seqfish, _load_seqfishplus
from scvi.dataset._built_in_data._smfish import _load_smfish
from scvi.dataset._built_in_data._dataset10X import _load_dataset10X


def pbmc_dataset(
    save_path: str = "data/",
    run_setup_anndata: bool = True,
    remove_extracted_data: bool = True,
):
    """Loads pbmc dataset.

    We considered scRNA-seq data from two batches of peripheral blood mononuclear cells (PBMCs) from a healthy donor
    (4K PBMCs and 8K PBMCs). We derived quality control metrics using the cellrangerRkit R package (v. 1.1.0).
    Quality metrics were extracted from CellRanger throughout the molecule specific information file. After filtering,
    we extract 12,039 cells with 10,310 sampled genes and get biologically meaningful clusters with the
    software Seurat. We then filter genes that we could not match with the bulk data used for differential
    expression to be left with g = 3346.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning
    remove_extracted_data
        If true, will remove the folder the data was extracted to

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = pbmc_dataset()
    """

    return _load_pbmc_dataset(
        save_path=save_path,
        run_setup_anndata=run_setup_anndata,
        remove_extracted_data=remove_extracted_data,
    )


def dataset10X(
    dataset_name: Optional[str] = None,
    filename: Optional[str] = None,
    save_path: str = "data/10X",
    url: str = None,
    is_filtered: bool = True,
    remove_extracted_data: bool = False,
    **scanpy_read_10x_kwargs,
):
    """Loads a file from `10x <http://cf.10xgenomics.com/>`_ website.
    Parameters
    ----------
    dataset_name
        Name of the dataset file. Has to be one of:
        "frozen_pbmc_donor_a", "frozen_pbmc_donor_b", "frozen_pbmc_donor_c", "fresh_68k_pbmc_donor_a",
        "cd14_monocytes", "b_cells", "cd34", "cd56_nk", "cd4_t_helper", "regulatory_t", "naive_t",
        "memory_t", "cytotoxic_t", "naive_cytotoxic", "pbmc8k", "pbmc4k", "t_3k", "t_4k", "neuron_9k",
        "pbmc_1k_protein_v3", "pbmc_10k_protein_v3", "malt_10k_protein_v3", "pbmc_1k_v2", "pbmc_1k_v3",
        "pbmc_10k_v3", "hgmm_1k_v2", "hgmm_1k_v3", "hgmm_5k_v3", "hgmm_10k_v3", "neuron_1k_v2",
        "neuron_1k_v3", "neuron_10k_v3", "heart_1k_v2", "heart_1k_v3", "heart_10k_v3".
    filename
        manual override of the filename to write to.
    save_path
        Location to use when saving/loading the data.
    url
        manual override of the download remote location.
        Note that we already provide urls for most 10X datasets,
        which are automatically formed only using the ``dataset_name``.
    type
        Either `filtered` data or `raw` data.
    remove_extracted_data
        Whether to remove extracted archives in the case of `.tar.gz` downloads.

    Examples
    --------
    >>> neuron = dataset10X("neuron_9k")
    """
    return _load_dataset10X(
        dataset_name=dataset_name,
        filename=filename,
        save_path=save_path,
        url=url,
        is_filtered=is_filtered,
        remove_extracted_data=remove_extracted_data,
        **scanpy_read_10x_kwargs,
    )


def smfish(
    save_path: str = "data/",
    use_high_level_cluster: bool = True,
    run_setup_anndata: bool = True,
):
    """Loads osmFISH data of mouse cortex cells from the Linarsson lab.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    use_high_level_cluster
        If True, use higher-level agglomerate clusters.
        The resulting cell types are "Astrocytes", "Endothelials", "Inhibitory",
        "Microglias", "Oligodendrocytes" and "Pyramidals".
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = smfish()
    """
    return _load_smfish(
        save_path=save_path,
        use_high_level_cluster=use_high_level_cluster,
        run_setup_anndata=run_setup_anndata,
    )


def seqfishplus(
    save_path: str = "data/",
    tissue_region="subventricular cortex",
    run_setup_anndata: bool = True,
):
    """seqFISH+ can image mRNAs for 10,000 genes in single cells—with high accuracy and
    sub-diffraction-limit resolution—in the cortex, subventricular zone
    and olfactory bulb of mouse brain

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    tissue_region
        Region of the mouse brain, Either "subventricular cortex" or "olfactory bulb"
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = seqfishplus()
    """
    return _load_seqfishplus(
        save_path=save_path,
        tissue_region=tissue_region,
        run_setup_anndata=run_setup_anndata,
    )


def seqfish(save_path: str = "data/", run_setup_anndata: bool = True):
    """seqfish dataset

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = seqfish()
    """
    return _load_seqfish(save_path=save_path, run_setup_anndata=run_setup_anndata)


def purified_pbmc_dataset(
    save_path: str = "data/",
    subset_datasets: Optional[List[str]] = None,
    run_setup_anndata: bool = True,
):
    """Purified PBMC dataset from: "Massively parallel digital transcriptional profiling of single cells".

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    subset_datasets
        index for subsetting the follwing list of datasets
        which are used to form the ``PurifiedPBMCDataset``:
        "cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic",
        "b_cells", "cd4_t_helper", "cd34", "cd56_nk", "cd14_monocytes".
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = purified_pbmc_dataset()
    """
    return _load_purified_pbmc_dataset(
        save_path=save_path,
        subset_datasets=subset_datasets,
        run_setup_anndata=run_setup_anndata,
    )


def prefrontalcortex_starmap(save_path: str = "data/", run_setup_anndata: bool = True):
    """Loads a starMAP dataset of 3,704 cells and 166 genes from
    the mouse pre-frontal cortex (Wang et al., 2018)

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = prefrontalcortex_starmap()
    """
    return _load_prefrontalcortex_starmap(save_path=save_path)


def frontalcortex_dropseq(save_path: str = "data/", run_setup_anndata: bool = True):
    """"Load the cells from the mouse frontal cortex sequenced by the Dropseq technology (Saunders et al., 2018)

    Load the 71639 annotated cells located in the frontal cortex of adult mouses among the 690,000 cells
    studied by (Saunders et al., 2018) using the Drop-seq method. We have a 71639*7611 gene expression matrix
    Among the 7611 genes, we offer the user to provide a list of genes to subsample from. If not provided,
    all genes are kept.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = frontalcortex_dropseq()
    """
    return _load_frontalcortex_dropseq(save_path=save_path)


def annotation_simulation(
    name: str, save_path: str = "data/", run_setup_anndata: bool = True
):
    """\
    Simulated datasets for scANVI tutorials

    Parameters
    ----------
    name
        One of "1", "2", or "3"
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = annontation_simulation("1")
    """
    return _load_annotation_simulation(
        name=name, save_path=save_path, run_setup_anndata=run_setup_anndata
    )


def retina(save_path: str = "data/", run_setup_anndata: bool = True):
    """Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells and
    13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from the author.
    We also extract their normalized data with Combat and use it for benchmarking.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = retina()
    """
    return _load_retina(save_path=save_path)


def mouse_ob_dataset(save_path: str = "data/", run_setup_anndata: bool = True):
    """Loads mouse ob dataset.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object with `.obs['batch']` and `.obs['labels']`

    Examples
    --------
    >>> adata = mouse_ob_dataset()
    """
    return _load_mouse_ob_dataset(
        save_path=save_path, run_setup_anndata=run_setup_anndata
    )


def breast_cancer_dataset(save_path: str = "data/", run_setup_anndata: bool = True):
    """Loads breast cancer dataset.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` object with `.obs['batch']` and `.obs['labels']`

    Examples
    --------
    >>> adata = breast_cancer_dataset()
    """
    return _load_breast_cancer_dataset(
        save_path=save_path, run_setup_anndata=run_setup_anndata
    )


def pbmcs_10x_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    run_setup_anndata: bool = True,
):
    """Filtered PBMCs from 10x Genomics profiled with RNA and protein

    Datasets were filtered for doublets and other outliers as in
    https://github.com/YosefLab/totalVI_reproducibility/blob/master/data/data_filtering_scripts/pbmc_10k/pbmc_10k.py

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` with `.obsm["protein_expression"]`

    Missing protein values are zero, and are identified during `AnnData` setup.

    Examples
    --------
    >>> adata = pbmcs_10x_cite_seq()
    """
    return _load_pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join=protein_join,
        run_setup_anndata=run_setup_anndata,
    )


def brainlarge_dataset(
    save_path: str = "data/",
    run_setup_anndata: bool = True,
    sample_size_gene_var: int = 10000,
    max_cells_to_keep: Optional[int] = None,
    n_genes_to_keep: int = 720,
    loading_batch_size: int = 100000,
):
    """Loads brain-large dataset.

    This dataset contains 1.3 million brain cells from
    `10x Genomics <https://support.10xgenomics.com/single-cell-gene-expression/datasets>`_.
    We randomly shuffle the data to get a 1M subset of cells and order genes by variance to retain first 10,000 and then 720 sampled variable genes.
    This dataset is then sampled multiple times in cells for the runtime and goodness-of-fit analysis.
    We report imputation scores on the 10k cells and 720 genes samples only.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning
    sample_size_gene_var
        Number of cells to use to estimate gene variances.
    max_cells_to_keep
        Maximum number of cells to keep.
    n_genes_to_keep
        Number of genes to keep, ordered by decreasing variance.
    loading_batch_size
        Number of cells to use for each chunk loaded.

    Examples
    --------
    >>> adata = brainlarge_dataset()

    Returns
    -------
    `AnnData` object

    """

    return _load_brainlarge_dataset(
        save_path=save_path,
        run_setup_anndata=run_setup_anndata,
        sample_size_gene_var=sample_size_gene_var,
        max_cells_to_keep=max_cells_to_keep,
        n_genes_to_keep=n_genes_to_keep,
        loading_batch_size=loading_batch_size,
    )


def cortex(save_path: str = "data/", run_setup_anndata: bool = True):
    """
    Loads cortex dataset.

    The
    `Mouse Cortex Cells dataset <https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt>`_
    contains 3005 mouse cortex cells and gold-standard labels for seven distinct cell types. Each cell type corresponds
    to a cluster to recover.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Examples
    --------
    >>> adata = cortex()

    Returns
    -------
    `AnnData` object
    """
    return _load_cortex(save_path, run_setup_anndata)


def synthetic_iid(
    batch_size: Optional[int] = 200,
    n_genes: Optional[int] = 100,
    n_proteins: Optional[int] = 100,
    n_batches: Optional[int] = 2,
    n_labels: Optional[int] = 3,
) -> anndata.AnnData:
    """Synthetic dataset with ZINB distributed RNA and NB distributed protein

    Each value is independently and identically distributed.

    Parameters
    ----------
    batch_size
        Number of cells per batch
    n_genes
        Number of genes
    n_proteins
        Number of proteins
    n_batches
        Number of batches
    n_labels
        Number of cell types

    Returns
    -------
    `AnnData` object

    Examples
    --------
    >>> adata = synthetic_iid()
    """

    return _generate_synthetic(
        batch_size=batch_size,
        n_genes=n_genes,
        n_proteins=n_proteins,
        n_batches=n_batches,
        n_labels=n_labels,
    )
