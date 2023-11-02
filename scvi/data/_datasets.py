from __future__ import annotations

import anndata

from scvi._types import AnnOrMuData

from ._built_in_data._brain_large import _load_brainlarge_dataset
from ._built_in_data._cellxgene import _load_cellxgene_dataset
from ._built_in_data._cite_seq import (
    _load_pbmc_seurat_v4_cite_seq,
    _load_pbmcs_10x_cite_seq,
    _load_spleen_lymph_cite_seq,
)
from ._built_in_data._cortex import _load_cortex
from ._built_in_data._csv import _load_breast_cancer_dataset, _load_mouse_ob_dataset
from ._built_in_data._dataset_10x import _load_dataset_10x
from ._built_in_data._heartcellatlas import _load_heart_cell_atlas_subsampled
from ._built_in_data._loom import (
    _load_annotation_simulation,
    _load_frontalcortex_dropseq,
    _load_prefrontalcortex_starmap,
    _load_retina,
)
from ._built_in_data._pbmc import _load_pbmc_dataset, _load_purified_pbmc_dataset
from ._built_in_data._smfish import _load_smfish
from ._built_in_data._synthetic import _generate_synthetic


def pbmc_dataset(
    save_path: str = "data/",
    remove_extracted_data: bool = True,
) -> anndata.AnnData:
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
    remove_extracted_data
        If true, will remove the folder the data was extracted to

    Returns
    -------
    AnnData with batch info (``.obs['batch']``), label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.pbmc_dataset()
    """
    return _load_pbmc_dataset(
        save_path=save_path,
        remove_extracted_data=remove_extracted_data,
    )


def dataset_10x(
    dataset_name: str | None = None,
    filename: str | None = None,
    save_path: str = "data/10X",
    url: str = None,
    return_filtered: bool = True,
    remove_extracted_data: bool = False,
    **scanpy_read_10x_kwargs,
) -> anndata.AnnData:
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
        "neuron_1k_v3", "neuron_10k_v3", "heart_1k_v2", "heart_1k_v3", "heart_10k_v3", 5k_pbmc_protein_v3",
        "5k_pbmc_protein_v3_nextgem", 1M_neurons".
    filename
        manual override of the filename to write to.
    save_path
        Location to use when saving/loading the data.
    url
        manual override of the download remote location.
        Note that we already provide urls for most 10X datasets,
        which are automatically formed only using the ``dataset_name``.
    return_filtered
        Either `filtered` data or `raw` data.
    remove_extracted_data
        Whether to remove extracted archives in the case of `.tar.gz` downloads.
    **scanpy_read_10x_kwargs
        Kwargs for scanpy's read_10x function

    Returns
    -------
    adata initialized with 10x data

    Examples
    --------
    >>> import scvi
    >>> neuron = scvi.data.dataset10X("neuron_9k")
    """
    return _load_dataset_10x(
        dataset_name=dataset_name,
        filename=filename,
        save_path=save_path,
        url=url,
        return_filtered=return_filtered,
        remove_extracted_data=remove_extracted_data,
        **scanpy_read_10x_kwargs,
    )


def cellxgene(
    url: str,
    filename: str | None = None,
    save_path: str = "data/",
    return_path: bool = False,
) -> anndata.AnnData | str:
    """Loads a file from `cellxgene <https://cellxgene.cziscience.com/>`_ portal.

    Parameters
    ----------
    url
        URL to cellxgene session
    filename
        manual override of the filename to write to.
    save_path
        Location to use when saving/loading the data.
    return_path
        Whether to return the path to the downloaded file.

    Returns
    -------
    adata initialized with cellxgene data
    """
    return _load_cellxgene_dataset(
        url=url,
        filename=filename,
        save_path=save_path,
        return_path=return_path,
    )


def smfish(
    save_path: str = "data/",
    use_high_level_cluster: bool = True,
) -> anndata.AnnData:
    """Loads osmFISH data of mouse cortex cells from the Linarsson lab.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    use_high_level_cluster
        If True, use higher-level agglomerate clusters.
        The resulting cell types are "Astrocytes", "Endothelials", "Inhibitory",
        "Microglias", "Oligodendrocytes" and "Pyramidals".

    Returns
    -------
    AnnData with batch info (``.obs['batch']``), label info (``.obs['labels']``),
    spatial info (``.obs['x_coord']``, ``.obs['y_coord']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.smfish()
    """
    return _load_smfish(
        save_path=save_path,
        use_high_level_cluster=use_high_level_cluster,
    )


def purified_pbmc_dataset(
    save_path: str = "data/",
    subset_datasets: list[str] | None = None,
) -> anndata.AnnData:
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

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.purified_pbmc_dataset()
    """
    return _load_purified_pbmc_dataset(
        save_path=save_path,
        subset_datasets=subset_datasets,
    )


def prefrontalcortex_starmap(save_path: str = "data/") -> anndata.AnnData:
    """Loads a starMAP dataset of mouse pre-frontal cortex (Wang et al., 2018).

    3,704 cells and 166 genes.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``), label info (``.obs['labels']``),
    spatial info (``.obs['x_coord']``, ``.obs['y_coord']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.prefrontalcortex_starmap()

    """
    return _load_prefrontalcortex_starmap(save_path=save_path)


def frontalcortex_dropseq(save_path: str = "data/") -> anndata.AnnData:
    """Load the cells from the mouse frontal cortex sequenced by the Dropseq technology (Saunders et al., 2018).

    Load the 71639 annotated cells located in the frontal cortex of adult mouses among the 690,000 cells
    studied by (Saunders et al., 2018) using the Drop-seq method. We have a 71639*7611 gene expression matrix
    Among the 7611 genes, we offer the user to provide a list of genes to subsample from. If not provided,
    all genes are kept.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.frontalcortex_dropseq()
    """
    return _load_frontalcortex_dropseq(save_path=save_path)


def annotation_simulation(name: str, save_path: str = "data/") -> anndata.AnnData:
    """Simulated datasets for scANVI tutorials.

    Parameters
    ----------
    name
        One of "1", "2", or "3"
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.annontation_simulation("1")

    """
    return _load_annotation_simulation(name=name, save_path=save_path)


def retina(save_path: str = "data/") -> anndata.AnnData:
    """Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells and
    13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from the author.
    We also extract their normalized data with Combat and use it for benchmarking.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> adata = retina()
    """
    return _load_retina(save_path=save_path)


def mouse_ob_dataset(save_path: str = "data/") -> anndata.AnnData:
    """Loads mouse ob dataset.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.mouse_ob_dataset()
    """
    return _load_mouse_ob_dataset(save_path=save_path)


def breast_cancer_dataset(save_path: str = "data/") -> anndata.AnnData:
    """Loads breast cancer dataset.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.breast_cancer_dataset()
    """
    return _load_breast_cancer_dataset(save_path=save_path)


def pbmcs_10x_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
) -> anndata.AnnData:
    """Filtered PBMCs from 10x Genomics profiled with RNA and protein.

    Datasets were filtered for doublets and other outliers as in
    https://github.com/YosefLab/totalVI_reproducibility/blob/master/data/data_filtering_scripts/pbmc_10k/pbmc_10k.py

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins

    Returns
    -------
    AnnData with batch info (``.obs['batch']``),
    and protein expression (``.obsm["protein_expression"]``)

    Missing protein values are zero, when ``protein_join == "outer`` and are identified during ``AnnData`` setup.

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.pbmcs_10x_cite_seq()
    """
    return _load_pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join=protein_join,
    )


def spleen_lymph_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    remove_outliers: bool = True,
) -> anndata.AnnData:
    """Immune cells from the murine spleen and lymph nodes :cite:p:`GayosoSteier21`.

    This dataset was used throughout the totalVI manuscript, and named SLN-all.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins
    remove_outliers
        Whether to remove clusters annotated as doublet or low quality

    Returns
    -------
    AnnData with batch info (``.obs['batch']``), label info (``.obs['cell_types']``),
    protein expression (``.obsm["protein_expression"]``), and tissue (``.obs['tissue']``).

    Missing protein values are zero, when ``protein_join == "outer`` and are identified during ``AnnData`` setup.

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.spleen_lymph_cite_seq()
    """
    return _load_spleen_lymph_cite_seq(
        save_path=save_path,
        protein_join=protein_join,
        remove_outliers=remove_outliers,
    )


def brainlarge_dataset(
    save_path: str = "data/",
    sample_size_gene_var: int = 10000,
    max_cells_to_keep: int | None = None,
    n_genes_to_keep: int = 720,
    loading_batch_size: int = 100000,
) -> anndata.AnnData:
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
    sample_size_gene_var
        Number of cells to use to estimate gene variances.
    max_cells_to_keep
        Maximum number of cells to keep.
    n_genes_to_keep
        Number of genes to keep, ordered by decreasing variance.
    loading_batch_size
        Number of cells to use for each chunk loaded.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.brainlarge_dataset()
    """
    return _load_brainlarge_dataset(
        save_path=save_path,
        sample_size_gene_var=sample_size_gene_var,
        max_cells_to_keep=max_cells_to_keep,
        n_genes_to_keep=n_genes_to_keep,
        loading_batch_size=loading_batch_size,
    )


def cortex(save_path: str = "data/") -> anndata.AnnData:
    """Loads cortex dataset.

    The
    `Mouse Cortex Cells dataset <https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt>`_
    contains 3005 mouse cortex cells and gold-standard labels for seven distinct cell types. Each cell type corresponds
    to a cluster to recover.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    AnnData with batch info (``.obs['batch']``) and label info (``.obs['labels']``)

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.cortex()
    """
    return _load_cortex(save_path)


def synthetic_iid(
    batch_size: int = 200,
    n_genes: int = 100,
    n_proteins: int = 100,
    n_regions: int = 100,
    n_batches: int = 2,
    n_labels: int = 3,
    dropout_ratio: float = 0.7,
    sparse_format: str | None = None,
    return_mudata: bool = False,
) -> AnnOrMuData:
    """Synthetic multimodal dataset.

    RNA and accessibility data are generated from a zero-inflated negative binomial,
    while protein data is generated from a negative binomial distribution. This dataset
    is just for testing purposes and not meant for modeling or research. Each value is
    independently and identically distributed.

    Parameters
    ----------
    batch_size
        The number of cells per batch such that the total number of cells in the data is
        `batch_size * n_batches`.
    n_genes
        The number of genes to generate.
    n_proteins
        The number of proteins to generate.
    n_regions
        The number of accessibility regions to generate.
    n_batches
        The number of batches to generate.
    n_labels
        The number of cell type labels, distributed uniformly across batches.
    sparse
        Whether to store ZINB generated data as a :class:`scipy.sparse.csr_matrix`.
    dropout_ratio
        The expected percentage of zeros artificially added into the data for RNA
        and accessibility data.
    sparse_format
        Whether to store RNA, accessibility, and protein data as sparse arrays. One of
        the following:

        * `None`: Store as a dense :class:`numpy.ndarray`.
        * `"csr_matrix"`: Store as a :class:`scipy.sparse.csr_matrix`.
        * `"csc_matrix"`: Store as a :class:`scipy.sparse.csc_matrix`.
    return_mudata
        Returns a :class:`~mudata.MuData` if `True`, else :class:`~anndata.AnnData`.

    Returns
    -------
    :class:`~anndata.AnnData` (if `return_mudata=False`) with the following fields:

    * `.obs["batch"]`: Categorical batch labels in the format `batch_{i}`.
    * `.obs["labels"]`: Categorical cell type labels in the format `label_{i}`.
    * `.obsm["protein_expression"]`: Protein expression matrix.
    * `.uns["protein_names"]`: Array of protein names.
    * `.obsm["accessibility"]`: Accessibility expression matrix.

    :class:`~mudata.MuData` (if `return_mudata=True`) with the following fields:

    * `.obs["batch"]`: Categorical batch labels in the format `batch_{i}`.
    * `.obs["labels"]`: Categorical cell type labels in the format `label_{i}`.
    * `.mod["rna"]`: RNA expression data.
    * `.mod["protein_expression"]`: Protein expression data.
    * `.mod["accessibility"]`: Accessibility expression data.

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.synthetic_iid()
    """
    if n_batches < 1:
        raise ValueError("`n_batches` must be greater than 0")
    if n_genes < 1:
        raise ValueError("`n_genes` must be greater than 0")

    return _generate_synthetic(
        batch_size=batch_size,
        n_genes=n_genes,
        n_proteins=n_proteins,
        n_regions=n_regions,
        n_batches=n_batches,
        n_labels=n_labels,
        dropout_ratio=dropout_ratio,
        sparse_format=sparse_format,
        return_mudata=return_mudata,
    )


def heart_cell_atlas_subsampled(
    save_path: str = "data/",
    remove_nuisance_clusters: bool = True,
) -> anndata.AnnData:
    """Combined single cell and single nuclei RNA-Seq data of 485K cardiac cells with annotations.

    Dataset was filtered down randomly to 20k cells using :meth:`~scanpy.pp.subsample`. The original
    data can be downloaded from https://www.heartcellatlas.org/#DataSources.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    remove_nuisance_clusters
        Remove doublets and unsassigned cells

    Returns
    -------
    AnnData

    Notes
    -----
    The data were filtered using the following sequence::

        >>> adata = anndata.read_h5ad(path_to_anndata)
        >>> bdata = sc.pp.subsample(adata, n_obs=20000, copy=True)
        >>> sc.pp.filter_genes(bdata, min_counts=3)
        >>> bdata.write_h5ad(path, compression="gzip")

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.heart_cell_atlas_subsampled()
    """
    return _load_heart_cell_atlas_subsampled(
        save_path=save_path,
        remove_nuisance_clusters=remove_nuisance_clusters,
    )


def pbmc_seurat_v4_cite_seq(
    save_path: str = "data/",
    apply_filters: bool = True,
    aggregate_proteins: bool = True,
    mask_protein_batches: int = 0,
) -> anndata.AnnData:
    """Dataset of PBMCs measured with CITE-seq (161764 cells).

    This dataset was first presented in the Seurat v4 paper:

    https://doi.org/10.1016/j.cell.2021.04.048

    It contains 8 volunteers in an HIV vaccine trial measured
    at 3 time points; thus, there are 24 batches in this dataset.


    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    apply_filters
        Apply filters at cell and protein level. At the cell level,
        this filters on protein library size, number proteins detected,
        percent mito, and removes cells labeled as doublets.
    aggregate_proteins
        Antibodies targeting the same surface protein are added together,
        and isotype controls are removed. See the source code for full details.
    mask_protein_batches
        Set proteins in this many batches to be all zero (considered missing
        for :class:`~scvi.model.TOTALVI`.). This improves transfer learning
        with this dataset.

    Returns
    -------
    AnnData

    Notes
    -----
    This is not the same exact dataset as can be downloaded from:

    https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

    This is due to the fact that the object linked in the tutorial above does
    not contain the actual UMI count data for RNA. UMI counts had to be separately
    downloaded from GEO (GSE164378). The counts in that object are an output of the
    scTransform method and should not be treated like UMI counts.

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.pbmc_seurat_v4_cite_seq()
    """
    return _load_pbmc_seurat_v4_cite_seq(
        save_path=save_path,
        apply_filters=apply_filters,
        aggregate_proteins=aggregate_proteins,
        mask_protein_batches=mask_protein_batches,
    )
