import numpy as np
from .dataset10X import Dataset10X
from .dataset import GeneExpressionDataset


class PbmcDataset(GeneExpressionDataset):
    r""" Loads pbmc dataset.

    We considered scRNA-seq data from two batches of peripheral blood mononuclear cells (PBMCs) from a healthy donor
    (4K PBMCs and 8K PBMCs). We derived quality control metrics using the cellrangerRkit R package (v. 1.1.0).
    Quality metrics were extracted from CellRanger throughout the molecule specific information file. After filtering,
    we extract 12,039 cells with 10,310 sampled genes and get biologically meaningful clusters with the
    software Seurat. We then filter genes that we could not match with the bulk data used for differential
    expression to be left with g = 3346.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = PbmcDataset()

    """
    def __init__(self, save_path='data/'):
        pbmc = GeneExpressionDataset.concat_datasets(
                    Dataset10X("pbmc8k", save_path=save_path),
                    Dataset10X("pbmc4k", save_path=save_path)
                )
        super(PbmcDataset, self).__init__(pbmc.X, pbmc.local_means, pbmc.local_vars,
                                          pbmc.batch_indices, pbmc.labels, pbmc.gene_names)


class PurifiedPBMCDataset(GeneExpressionDataset):
    r""" The purified PBMC dataset from: "Massively parallel digital transcriptional profiling of single cells".

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = PurifiedPBMCDataset()

    """

    def __init__(self, save_path='data/', filter_cell_types=None):
        cell_types = np.array(["cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic",
                               "b_cells", "cd4_t_helper", "cd34", "cd56_nk", "cd14_monocytes"])
        if filter_cell_types:  # filter = np.arange(6) - for T cells:  np.arange(4) for T/CD4 cells
            cell_types = cell_types[np.array(filter_cell_types)]

        datasets = []
        for cell_type in cell_types:
            dataset = Dataset10X(cell_type, save_path=save_path)
            dataset.cell_types = np.array([cell_type])
            datasets += [dataset]

        pbmc = GeneExpressionDataset.concat_datasets(*datasets, shared_batches=True)
        pbmc.subsample_genes(subset_genes=(np.array(pbmc.X.sum(axis=0)) > 0).ravel())
        super(PurifiedPBMCDataset, self).__init__(pbmc.X, pbmc.local_means, pbmc.local_vars,
                                                  pbmc.batch_indices, pbmc.labels,
                                                  gene_names=pbmc.gene_names, cell_types=pbmc.cell_types)
