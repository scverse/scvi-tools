from unittest import TestCase

import scvi
import tarfile
import os
import scanpy as sc
from .utils import unsupervised_training_one_epoch


class TestDataset10X(TestCase):
    def test_populate_and_train_one_v1(self):
        dataset = scvi.dataset.dataset10X(
            dataset_name="cd4_t_helper",
            remove_extracted_data=True,
            save_path="tests/data/10X",
        )
        scvi.dataset.setup_anndata(dataset)
        unsupervised_training_one_epoch(dataset)

    def test_brain_small(self):
        dataset = scvi.dataset.dataset10X(
            dataset_name="neuron_9k",
            save_path="tests/data/10X",
            remove_extracted_data=True,
        )
        scvi.dataset.setup_anndata(dataset)
        unsupervised_training_one_epoch(dataset)

    def test_pbmc_cite(self):
        file_path = (
            "tests/data/10X/pbmc_10k_protein_v3/filtered_feature_bc_matrix.tar.gz"
        )
        save_path = "tests/data/10X/pbmc_10k_protein_v3/"
        tar = tarfile.open(file_path, "r:gz")
        tar.extractall(path=save_path)
        tar.close()
        dataset = sc.read_10x_mtx(
            os.path.join(save_path, "filtered_feature_bc_matrix"), gex_only=False
        )
        scvi.dataset.organize_cite_seq_10x(dataset)
        scvi.dataset.setup_anndata(
            dataset, protein_expression_obsm_key="protein_expression"
        )
        unsupervised_training_one_epoch(dataset)
