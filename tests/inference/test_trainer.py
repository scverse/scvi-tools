from unittest import TestCase
import numpy as np
from scvi.dataset import GeneExpressionDataset, CellMeasurement
from scvi.inference import UnsupervisedTrainer, JVAETrainer, TotalTrainer
from scvi.models import VAE, JVAE, TOTALVI, Classifier


class TestTrainer(TestCase):
    def test_special_dataset_size(self):
        gene_dataset = GeneExpressionDataset()
        x = np.random.randint(1, 100, (17 * 2, 10))
        y = np.random.randint(1, 100, (17 * 2, 10))
        gene_dataset.populate_from_data(x)
        protein_data = CellMeasurement(
            name="protein_expression",
            data=y,
            columns_attr_name="protein_names",
            columns=np.arange(10),
        )
        gene_dataset.initialize_cell_measurement(protein_data)

        # Test UnsupervisedTrainer
        vae = VAE(
            gene_dataset.nb_genes,
            n_batch=gene_dataset.n_batches,
            n_labels=gene_dataset.n_labels,
        )
        trainer = UnsupervisedTrainer(
            vae,
            gene_dataset,
            train_size=0.5,
            use_cuda=False,
            data_loader_kwargs={"batch_size": 8},
        )
        trainer.train(n_epochs=1)

        # Test JVATrainer
        jvae = JVAE(
            [gene_dataset.nb_genes, gene_dataset.nb_genes],
            gene_dataset.nb_genes,
            [slice(None)] * 2,
            ["zinb", "zinb"],
            [True, True],
            n_batch=1,
        )
        cls = Classifier(gene_dataset.nb_genes, n_labels=2, logits=True)
        trainer = JVAETrainer(
            jvae,
            cls,
            [gene_dataset, gene_dataset],
            train_size=0.5,
            use_cuda=False,
            data_loader_kwargs={"batch_size": 8},
        )
        trainer.train(n_epochs=1)

        totalvae = TOTALVI(gene_dataset.nb_genes, len(gene_dataset.protein_names))
        trainer = TotalTrainer(
            totalvae,
            gene_dataset,
            train_size=0.5,
            use_cuda=False,
            data_loader_kwargs={"batch_size": 8},
        )
        trainer.train(n_epochs=1)
