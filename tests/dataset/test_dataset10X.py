from unittest import TestCase

import scvi
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
        dataset = scvi.dataset.dataset10X(
            dataset_name="pbmc_10k_v3",
            save_path="tests/data/10X",
            remove_extracted_data=True,
        )
        scvi.dataset.organize_cite_seq_10x(dataset)
        scvi.dataset.setup_anndata(
            dataset, protein_expression_obsm_key="protein_expression"
        )
        unsupervised_training_one_epoch(dataset)
