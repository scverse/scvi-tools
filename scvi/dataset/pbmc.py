import os
import pickle
from typing import Union, List

import numpy as np
import pandas as pd

from scvi.dataset.dataset import DownloadableDataset, remap_categories
from scvi.dataset.dataset10X import Dataset10X


class PbmcDataset(DownloadableDataset):
    """Loads pbmc dataset.

    We considered scRNA-seq data from two batches of peripheral blood mononuclear cells (PBMCs) from a healthy donor
    (4K PBMCs and 8K PBMCs). We derived quality control metrics using the cellrangerRkit R package (v. 1.1.0).
    Quality metrics were extracted from CellRanger throughout the molecule specific information file. After filtering,
    we extract 12,039 cells with 10,310 sampled genes and get biologically meaningful clusters with the
    software Seurat. We then filter genes that we could not match with the bulk data used for differential
    expression to be left with g = 3346.

    :param save_path: Location to use when saving/loading the Pbmc metadata.
    :param save_path_10X: Location to use when saving/loading the underlying 10X datasets.
    :param remove_extracted_data: Whether to remove extracted archives after populating the dataset.
    :param delayed_populating: Switch for delayed populating mechanism.

    Examples:
        >>> gene_dataset = PbmcDataset()
    """

    def __init__(
        self,
        save_path: str = "data/",
        save_path_10X: str = None,
        remove_extracted_data: bool = False,
        delayed_populating: bool = False,
    ):
        self.save_path_10X = save_path_10X if save_path_10X is not None else save_path
        self.remove_extracted_data = remove_extracted_data
        self.barcodes = None
        super().__init__(
            urls=[
                "https://github.com/YosefLab/scVI-data/raw/master/gene_info.csv",
                "https://github.com/YosefLab/scVI-data/raw/master/pbmc_metadata.pickle",
            ],
            filenames=["gene_info_pbmc.csv", "pbmc_metadata.pickle"],
            save_path=save_path,
            delayed_populating=delayed_populating,
        )
        # this downloads the necessary file for a future call to populate
        if delayed_populating:
            Dataset10X("pbmc8k", save_path=self.save_path_10X, delayed_populating=True)
            Dataset10X("pbmc4k", save_path=self.save_path_10X, delayed_populating=True)

    def populate(self):
        self.de_metadata = pd.read_csv(
            os.path.join(self.save_path, "gene_info_pbmc.csv"), sep=","
        )
        pbmc_metadata = pickle.load(
            open(os.path.join(self.save_path, "pbmc_metadata.pickle"), "rb")
        )
        datasets = [
            Dataset10X(
                "pbmc8k",
                save_path=self.save_path_10X,
                remove_extracted_data=self.remove_extracted_data,
            ),
            Dataset10X(
                "pbmc4k",
                save_path=self.save_path_10X,
                remove_extracted_data=self.remove_extracted_data,
            ),
        ]
        self.populate_from_datasets(datasets)
        # filter cells according to barcodes
        dict_barcodes = dict(zip(self.barcodes, np.arange(len(self.barcodes))))
        subset_cells = []
        barcodes_metadata = (
            pbmc_metadata["barcodes"].index.values.ravel().astype(np.str)
        )
        for barcode in barcodes_metadata:
            if (
                barcode in dict_barcodes
            ):  # barcodes with end -11 filtered on 10X website (49 cells)
                subset_cells += [dict_barcodes[barcode]]
        self.update_cells(subset_cells=np.asarray(subset_cells))
        idx_metadata = np.asarray(
            [not barcode.endswith("11") for barcode in barcodes_metadata], dtype=np.bool
        )
        labels = pbmc_metadata["clusters"][idx_metadata].reshape(-1, 1)[: len(self)]
        self.labels, self.n_labels = remap_categories(labels)
        self.cell_types = pbmc_metadata["list_clusters"][: self.n_labels]

        genes_to_keep = list(
            self.de_metadata["ENSG"].values
        )  # only keep the genes for which we have de data
        difference = list(
            set(genes_to_keep).difference(set(self.gene_names))
        )  # Non empty only for unit tests
        for gene in difference:
            genes_to_keep.remove(gene)
        self.filter_genes_by_attribute(genes_to_keep)
        self.de_metadata = self.de_metadata.head(
            len(genes_to_keep)
        )  # this would only affect the unit tests
        self.design = pbmc_metadata["design"][idx_metadata]
        self.raw_qc = pbmc_metadata["raw_qc"][idx_metadata]
        self.qc_names = self.raw_qc.columns
        self.qc = self.raw_qc.values

        self.qc_pc = pbmc_metadata["qc_pc"][idx_metadata]
        self.normalized_qc = pbmc_metadata["normalized_qc"][idx_metadata]


class PurifiedPBMCDataset(DownloadableDataset):
    """Purified PBMC dataset from: "Massively parallel digital transcriptional profiling of single cells".

    :param subset_datasets: index for subsetting the follwing list of datasets
        which are used to form the ``PurifiedPBMCDataset``:
        "cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic",
        "b_cells", "cd4_t_helper", "cd34", "cd56_nk", "cd14_monocytes".

    Examples:
        >>> gene_dataset = PurifiedPBMCDataset()
    """

    def __init__(
        self,
        save_path: str = "data/",
        subset_datasets: Union[List[int], np.ndarray] = None,
        remove_extracted_data: bool = False,
        delayed_populating: bool = False,
    ):
        self.dataset_names = np.asarray(
            [
                "cd4_t_helper",
                "regulatory_t",
                "naive_t",
                "memory_t",
                "cytotoxic_t",
                "naive_cytotoxic",
                "b_cells",
                "cd4_t_helper",
                "cd34",
                "cd56_nk",
                "cd14_monocytes",
            ]
        )
        subset_datasets = subset_datasets if subset_datasets else slice(None)
        self.remove_extracted_data = remove_extracted_data
        self.dataset_names = self.dataset_names[subset_datasets]
        super().__init__(save_path=save_path, delayed_populating=delayed_populating)

        if delayed_populating:
            for dataset_name in self.dataset_names:
                Dataset10X(
                    dataset_name,
                    save_path=save_path,
                    delayed_populating=delayed_populating,
                )

        self.filter_genes_by_count()
        self.filter_cells_by_count()

    def populate(self):
        datasets = []
        for dataset_name in self.dataset_names:
            dataset = Dataset10X(
                dataset_name,
                save_path=self.save_path,
                remove_extracted_data=self.remove_extracted_data,
            )
            dataset.initialize_mapped_attribute(
                "labels", "cell_types", np.asarray([dataset_name], dtype=np.str)
            )
            datasets += [dataset]
        self.populate_from_datasets(datasets)
