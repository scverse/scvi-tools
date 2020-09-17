import scvi

from unittest import TestCase
from .utils import unsupervised_training_one_epoch


class TestPbmcDataset(TestCase):
    def test_populate(self):
        dataset = scvi.data.pbmc_dataset(
            save_path="tests/data/10X",
            remove_extracted_data=True,
            run_setup_anndata=True,
        )
        unsupervised_training_one_epoch(dataset)


class TestLoomDataset(TestCase):
    def test_retina_load_train_one(self):
        dataset = scvi.data.retina(save_path="tests/data")
        scvi.data.setup_anndata(dataset, batch_key="batch")
        unsupervised_training_one_epoch(dataset)

    def test_pfc_starmap_load_train_one(self):
        gene_dataset = scvi.data.prefrontalcortex_starmap(save_path="tests/data")
        scvi.data.setup_anndata(gene_dataset)
        unsupervised_training_one_epoch(gene_dataset)

    def test_fc_dropseq_load_train_one(self):
        gene_dataset = scvi.data.frontalcortex_dropseq(save_path="tests/data")
        scvi.data.setup_anndata(gene_dataset)
        unsupervised_training_one_epoch(gene_dataset)

    def test_smfish_load_train_one(self):
        gene_dataset = scvi.data.smfish(save_path="tests/data")
        scvi.data.setup_anndata(gene_dataset)
        unsupervised_training_one_epoch(gene_dataset)


class TestSeqfishDataset(TestCase):
    def test_populate(self):
        dataset = scvi.data.seqfish(save_path="tests/data")
        scvi.data.setup_anndata(dataset)
        unsupervised_training_one_epoch(dataset)


class TestSeqFishPlusDataset(TestCase):
    def test_populate(self):
        for tissue_region in ["subventricular cortex", "olfactory bulb"]:
            dataset = scvi.data.seqfishplus(
                tissue_region=tissue_region, save_path="tests/data"
            )
            scvi.data.setup_anndata(dataset)
            unsupervised_training_one_epoch(dataset)


class TestSyntheticDataset(TestCase):
    def test_iid(self):
        dataset = scvi.data.synthetic_iid(batch_size=10, n_genes=10)
        scvi.data.setup_anndata(dataset)
        unsupervised_training_one_epoch(dataset)


class TestCortexDataset(TestCase):
    def test_populate(self):
        adata = scvi.data.cortex(save_path="tests/data")
        scvi.data.setup_anndata(adata, labels_key="cell_type")
        unsupervised_training_one_epoch(adata)

    # def test_variance_and_order_and_size(self):
    #     to_keep = ["THY1", "sst", "Tomem2", "Crhbp"]
    #     total_genes = 10
    #     dataset_full = CortexDataset(save_path="tests/data", total_genes=None)
    #     dataset_small = CortexDataset(
    #         save_path="tests/data", genes_to_keep=to_keep, total_genes=total_genes
    #     )
    #     self.assertListEqual(dataset_small.gene_names[:4].tolist(), to_keep)

    #     small_variance = np.std(dataset_small.X[:, 4:], axis=0).argsort()[::-1]
    #     self.assertListEqual(small_variance.tolist(), list(range(6)))

    #     full_variance = np.std(dataset_full.X, axis=0).argsort()[::-1]
    #     variable_genes_all = dataset_full.gene_names[full_variance]
    #     genes_truth = (to_keep + [g for g in variable_genes_all if g not in to_keep])[
    #         :total_genes
    #     ]
    #     self.assertListEqual(dataset_small.gene_names.tolist(), genes_truth)


class TestBrainLargeDataset(TestCase):
    def test_populate(self):
        adata = scvi.data.brainlarge_dataset(
            save_path="tests/data",
            sample_size_gene_var=10,
            n_genes_to_keep=10,
            max_cells_to_keep=128,
        )
        unsupervised_training_one_epoch(adata)


class TestCsvDataset(TestCase):
    def test_breast_cancer(self):
        adata = scvi.data.breast_cancer_dataset(
            save_path="tests/data",
        )
        unsupervised_training_one_epoch(adata)

    def test_mouse_ob(self):
        adata = scvi.data.mouse_ob_dataset(
            save_path="tests/data",
        )
        unsupervised_training_one_epoch(adata)
