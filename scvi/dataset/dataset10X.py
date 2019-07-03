import logging
import os
import pickle
import tarfile
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.io as sp_io
import shutil
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import DownloadableDataset, _download

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
        If False, data is cast to sparse using ``scipy.sparse.csr_matrix``.
    :param gene_column: column in which to find gene names in the corresponding `.tsv` file.
    :param remove_extracted_data: Whether to remove extracted archives after populating the dataset.

    Examples:
        >>> tenX_dataset = Dataset10X("neuron_9k")

    .. _10x:
        http://cf.10xgenomics.com/
    """

    def __init__(
        self,
        dataset_name: str = None,
        filename: str = None,
        url: str = None,
        save_path: str = "data/10X",
        type: str = "filtered",
        dense: bool = False,
        gene_column: int = 0,
        remove_extracted_data: bool = False,
        delayed_populating: bool = False,
    ):
        self.barcodes = None
        self.dense = dense
        self.gene_column = gene_column
        self.remove_extracted_data = remove_extracted_data

        # form data url and filename unless manual override
        if dataset_name is not None:
            if url is not None:
                logger.warning("dataset_name provided, manual url is disregarded.")
            if filename is not None:
                logger.warning("dataset_name provided, manual filename is disregarded.")
            group = dataset_to_group[dataset_name]
            url_skeleton = group_to_url_skeleton[group]
            url = url_skeleton.format(group, dataset_name, dataset_name, type)
            filename = "%s_gene_bc_matrices.tar.gz" % type
            save_path = os.path.join(save_path, dataset_name)
        elif filename is not None and url is not None:
            logger.debug("Loading 10X dataset with custom url and filename")
        elif filename is not None and url is None:
            logger.debug("Loading local 10X dataset with custom filename")
        else:
            logger.debug("Loading extracted local 10X dataset with custom filename")
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing dataset")

        was_extracted = False
        if len(self.filenames) > 0:
            file_path = os.path.join(self.save_path, self.filenames[0])
            if not os.path.exists(file_path[:-7]):  # nothing extracted yet
                if tarfile.is_tarfile(file_path):
                    logger.info("Extracting tar file")
                    tar = tarfile.open(file_path, "r:gz")
                    tar.extractall(path=self.save_path)
                    was_extracted = True
                    tar.close()

        # get exact path, for robustness to changes is the 10X storage logic
        path_to_data, suffix = self.find_path_to_data()

        # get filenames, according to 10X storage logic
        if suffix == "":
            gene_filename = "genes.tsv"
        else:
            gene_filename = "features.tsv.gz"
        barcode_filename = "barcodes.tsv" + suffix

        genes_info = pd.read_csv(
            os.path.join(path_to_data, gene_filename), sep="\t", header=None
        )
        gene_names = genes_info.values[:, self.gene_column].astype(np.str)

        cell_attributes_dict = None
        batch_indices = None
        if os.path.exists(os.path.join(path_to_data, barcode_filename)):
            barcodes = pd.read_csv(
                os.path.join(path_to_data, barcode_filename), sep="\t", header=None
            )
            cell_attributes_dict = {
                "barcodes": np.squeeze(np.asarray(barcodes, dtype=str))
            }
            # As of 07/01, 10X barcodes have format "%s-%d" where the digit is a batch index starting at 1
            batch_indices = np.asarray(
                [barcode.split("-")[-1] for barcode in cell_attributes_dict["barcodes"]]
            )
            batch_indices = batch_indices.astype(np.int64) - 1

        matrix_filename = "matrix.mtx" + suffix
        expression_data = sp_io.mmread(os.path.join(path_to_data, matrix_filename)).T
        if self.dense:
            expression_data = expression_data.A
        else:
            expression_data = csr_matrix(expression_data)
        logger.info("Finished preprocessing dataset")

        self.populate_from_data(
            X=expression_data,
            batch_indices=batch_indices,
            gene_names=gene_names,
            cell_attributes_dict=cell_attributes_dict,
        )
        self.filter_cells_by_count()

        # cleanup if required
        if was_extracted and self.remove_extracted_data:
            logger.info("Removing extracted data at %s" % file_path[:-7])
            shutil.rmtree(file_path[:-7])

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
        raise FileNotFoundError(
            "No matrix.mtx(.gz) found in path (%s)." % self.save_path
        )


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

    def __init__(
        self,
        save_path: str = "data/",
        save_path_10X: str = None,
        delayed_populating: bool = False,
        remove_extracted_data: bool = False,
    ):
        super().__init__(
            dataset_name="neuron_9k",
            save_path=save_path_10X,
            remove_extracted_data=remove_extracted_data,
            delayed_populating=delayed_populating,
        )
        _download(
            filename="brain_small_metadata.pickle",
            url="https://github.com/YosefLab/scVI-data/raw/master/brain_small_metadata.pickle",
            save_path=save_path,
        )
        metadata = pickle.load(
            open(os.path.join(save_path, "brain_small_metadata.pickle"), "rb")
        )
        self.labels = metadata["clusters"].loc[self.barcodes] - 1
        raw_qc = metadata["raw_qc"].loc[self.barcodes]
        self.initialize_cell_attribute("raw_qc", raw_qc)
        self.initialize_cell_attribute("qc", raw_qc.values)
        self.qc_names = raw_qc.columns
