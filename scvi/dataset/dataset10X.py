import logging
import os
import pickle
import tarfile
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.io as sp_io
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)

available_datasets = {
    "1.1.0": [
        "frozen_pbmc_donor_a",
        "frozen_pbmc_donor_b",
        "frozen_pbmc_donor_c",
        "fresh_68k_pbmc_donor_a",
        "cd14_monocytes",
        "b_cells",
        "cd34",
        "cd56_nk",
        "cd4_t_helper",
        "regulatory_t",
        "naive_t",
        "memory_t",
        "cytotoxic_t",
        "naive_cytotoxic",
    ],
    "2.1.0": ["pbmc8k", "pbmc4k", "t_3k", "t_4k", "neuron_9k"],
    "3.0.0": [
        "pbmc_1k_protein_v3",
        "pbmc_10k_protein_v3",
        "malt_10k_protein_v3",
        "pbmc_1k_v2",
        "pbmc_1k_v3",
        "pbmc_10k_v3",
        "hgmm_1k_v2",
        "hgmm_1k_v3",
        "hgmm_5k_v3",
        "hgmm_10k_v3",
        "neuron_1k_v2",
        "neuron_1k_v3",
        "neuron_10k_v3",
        "heart_1k_v2",
        "heart_1k_v3",
        "heart_10k_v3",
    ],
}

dataset_to_group = dict(
    [
        (dataset_name, group)
        for group, list_datasets in available_datasets.items()
        for dataset_name in list_datasets
    ]
)
group_to_url_skeleton = {
    "1.1.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_gene_bc_matrices.tar.gz",
    "2.1.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_gene_bc_matrices.tar.gz",
    "3.0.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_feature_bc_matrix.tar.gz",
}
available_specification = ["filtered", "raw"]


class Dataset10X(DownloadableDataset):
    """Loads a file from `10x`_ website.

    :param dataset_name: Name of the dataset file. Has to be one of:
        "frozen_pbmc_donor_a", "frozen_pbmc_donor_b", "frozen_pbmc_donor_c", "fresh_68k_pbmc_donor_a",
        "cd14_monocytes", "b_cells", "cd34", "cd56_nk", "cd4_t_helper", "regulatory_t", "naive_t",
        "memory_t", "cytotoxic_t", "naive_cytotoxic", "pbmc8k", "pbmc4k", "t_3k", "t_4k", "neuron_9k",
        "pbmc_1k_protein_v3", "pbmc_10k_protein_v3", "malt_10k_protein_v3", "pbmc_1k_v2", "pbmc_1k_v3",
        "pbmc_10k_v3", "hgmm_1k_v2", "hgmm_1k_v3", "hgmm_5k_v3", "hgmm_10k_v3", "neuron_1k_v2",
        "neuron_1k_v3", "neuron_10k_v3", "heart_1k_v2", "heart_1k_v3", "heart_10k_v3".
    :param filename: manual override of the filename to write to.
    :param url: manual override of the download location.
        Note that we already provide urls for most 10X datasets,
        which are automatically formed only using the ``dataset_name``.
    :param save_path: Save path of the dataset.
    :param type: Either `filtered` data or `raw` data.
    :param dense: Whether to load as dense or sparse.
    :param gene_column: column in which to find gene names in the corresponding `.tsv` file.
    :param remote: Whether the 10X dataset is to be downloaded from the website
        or whether it is a local dataset, if remote is False then ``os.path.join(save_path, filename)``
        must be the path to the directory that contains matrix.mtx(.gz) and genes(features).tsv(.gz) files.

    Examples:
        >>> tenX_dataset = Dataset10X("neuron_9k")

    .. _10x:
        http://cf.10xgenomics.com/
    """

    def __init__(
        self,
        dataset_name: str = None,
        filename: str= None,
        url: str = None,
        save_path: str = "data/",
        type: str = "filtered",
        dense: bool = False,
        gene_column: int = 0,
        delayed_populating: bool = False,
    ):
        self.barcodes = None
        self.gene_column = gene_column
        self.dense = dense
        # form data url and filename unless manual override
        if url is None:
            group = dataset_to_group[dataset_name]
            url_skeleton = group_to_url_skeleton[group]
            url = url_skeleton.format(group, dataset_name, dataset_name, type)
            filename = "%s_gene_bc_matrices.tar.gz" % type

        super().__init__(
            urls=url,
            filenames=filename,
            save_path=os.path.join(save_path, "10X/%s/" % dataset_name),
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing dataset")
        file_path = os.path.join(self.save_path, self.filenames[0])
        if tarfile.is_tarfile(file_path):
            if not os.path.exists(file_path[:-6]):  # nothing extracted yet
                logger.info("Extracting tar file")
                tar = tarfile.open(file_path, "r:gz")
                tar.extractall(path=self.save_path)
                tar.close()
        path, suffix = self.find_path_to_data()
        # switch case corresponding to a change in 10X storing behaviour
        if suffix == "":
            gene_filename = "genes.tsv"
        else:
            gene_filename = "features.tsv.gz"
        genes_info = pd.read_csv(
            os.path.join(path, gene_filename), sep="\t", header=None
        )
        gene_names = genes_info.values[:, self.gene_column].astype(np.str)
        barcode_filename = "barcodes.tsv" + suffix
        cell_attributes_dict = None
        if os.path.exists(os.path.join(path, barcode_filename)):
            barcodes = pd.read_csv(
                os.path.join(path, barcode_filename), sep="\t", header=None
            )
            cell_attributes_dict = {"barcodes": barcodes}
        matrix_filename = "matrix.mtx" + suffix
        expression_data = sp_io.mmread(os.path.join(path, matrix_filename)).T
        if self.dense:
            expression_data = expression_data.A
        else:
            expression_data = csr_matrix(expression_data)
        logger.info("Finished preprocessing dataset")

        self.populate_from_data(
            X=expression_data,
            gene_names=gene_names,
            cell_attributes_dict=cell_attributes_dict,
        )
        self.filter_cells_by_count()

    def find_path_to_data(self) -> Tuple[str, str]:
        """Returns exact path for the data in the archive.

        This is required because 10X doesn't have a consistent way of storing their data.
        Additionally, the function returns whether the data is stored in compressed format.

        :return: path in which files are contains and their suffix if compressed.
        """
        for root, subdirs, files in os.walk(self.save_path):
            contains_mat = [
                filename == "matrix.mtx" or filename == "matrix.mtx.gz"
                for filename in files
            ]
            contains_mat = np.array(contains_mat).any()
            if contains_mat:
                is_tar = files[0][-3:] == ".gz"
                suffix = ".gz" if is_tar else ""
                return root, suffix
        raise FileNotFoundError("No matrix.mtx(.gz) found in path (%s)." % self.save_path)


class BrainSmallDataset(Dataset10X):
    """This dataset consists in 9,128 mouse brain cells profiled using `10x Genomics`.

    It is used as a complement of PBMC for our study of zero abundance
    and quality control metrics correlation with our generative posterior parameters.

    We derived quality control metrics using the cellrangerRkit R package (v.1.1.0).
    Quality metrics were extracted from CellRanger throughout the molecule specific information file.
    We kept the top 3000 genes by variance. We used the clusters provided by cellRanger
    for the correlation analysis of zero probabilities.

    Examples:
        >>> gene_dataset = BrainSmallDataset()

    .. _10x Genomics:
        https://support.10xgenomics.com/single-cell-gene-expression/datasets
    """
    def __init__(self, save_path: str = "data/", delayed_populating: bool = False):
        super().__init__(
            filename="brain_small_metadata.pickle",
            url="https://github.com/YosefLab/scVI-data/raw/master/brain_small_metadata.pickle",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )
        metadata = pickle.load(
            open(os.path.join(self.save_path, "brain_small_metadata.pickle"), "rb")
        )
        self.labels = metadata["clusters"].loc[self.barcodes.values.ravel()] - 1
        raw_qc = metadata["raw_qc"].loc[self.barcodes.values.ravel()]
        self.initialize_cell_attribute("raw_qc", raw_qc)
        self.initialize_cell_attribute("qc", raw_qc.values)
        self.qc_names = raw_qc.columns
