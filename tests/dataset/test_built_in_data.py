from unittest import TestCase

import pytest

import scvi

from .utils import unsupervised_training_one_epoch


class TestPbmcDataset(TestCase):
    def test_populate(self):
        dataset = scvi.data.pbmc_dataset(
            save_path="tests/data/10X",
            remove_extracted_data=True,
        )
        unsupervised_training_one_epoch(
            dataset,
        )


class TestLoomDataset(TestCase):
    def test_retina_load_train_one(self):
        dataset = scvi.data.retina(save_path="tests/data")
        unsupervised_training_one_epoch(dataset, batch_key="batch")

    def test_pfc_starmap_load_train_one(self):
        gene_dataset = scvi.data.prefrontalcortex_starmap(save_path="tests/data")
        unsupervised_training_one_epoch(gene_dataset)

    def test_fc_dropseq_load_train_one(self):
        gene_dataset = scvi.data.frontalcortex_dropseq(save_path="tests/data")
        unsupervised_training_one_epoch(gene_dataset)

    def test_smfish_load_train_one(self):
        gene_dataset = scvi.data.smfish(save_path="tests/data")
        unsupervised_training_one_epoch(gene_dataset)


class TestSeqfishDataset(TestCase):
    def test_populate(self):
        dataset = scvi.data.seqfish(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)


class TestSeqFishPlusDataset(TestCase):
    def test_populate(self):
        for tissue_region in ["subventricular cortex", "olfactory bulb"]:
            dataset = scvi.data.seqfishplus(
                tissue_region=tissue_region, save_path="tests/data"
            )
            unsupervised_training_one_epoch(dataset)


class TestSyntheticDataset(TestCase):
    def test_iid(self):
        dataset = scvi.data.synthetic_iid(batch_size=10, n_genes=10)
        unsupervised_training_one_epoch(dataset)


class TestCortexDataset(TestCase):
    def test_populate(self):
        adata = scvi.data.cortex(save_path="tests/data")
        unsupervised_training_one_epoch(adata, labels_key="cell_type")


class TestBrainLargeDataset(TestCase):
    def test_populate(self):
        adata = scvi.data.brainlarge_dataset(
            save_path="tests/data",
            sample_size_gene_var=10,
            n_genes_to_keep=10,
            max_cells_to_keep=128,
        )
        unsupervised_training_one_epoch(
            adata,
        )


class TestCsvDataset(TestCase):
    def test_breast_cancer(self):
        adata = scvi.data.breast_cancer_dataset(
            save_path="tests/data",
        )
        unsupervised_training_one_epoch(
            adata,
        )

    def test_mouse_ob(self):
        adata = scvi.data.mouse_ob_dataset(
            save_path="tests/data",
        )
        unsupervised_training_one_epoch(
            adata,
        )


@pytest.mark.internet
def test_download_spleen_lymph_data(save_path):
    scvi.data.spleen_lymph_cite_seq(save_path=save_path)
    scvi.data.spleen_lymph_cite_seq(save_path=save_path, protein_join="outer")


@pytest.mark.internet
def test_download_heart_cell_atlas(save_path):
    scvi.data.heart_cell_atlas_subsampled(save_path=save_path)


@pytest.mark.internet
def test_download_seurat_v4_pbmc(save_path):
    scvi.data.pbmc_seurat_v4_cite_seq(save_path=save_path, mask_protein_batches=5)
