#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np

from scvi.benchmark import all_benchmarks, benchmark, benchamrk_fish_scrna
from scvi.dataset import BrainLargeDataset, CortexDataset, RetinaDataset, BrainSmallDataset, HematoDataset, \
    LoomDataset, AnnDataset, CsvDataset, CiteSeqDataset, CbmcDataset, PbmcDataset, SyntheticDataset, \
    SeqfishDataset, SmfishDataset, BreastCancerDataset, MouseOBDataset, \
    GeneExpressionDataset, PurifiedPBMCDataset
from scvi.inference import JointSemiSupervisedTrainer, AlternateSemiSupervisedTrainer, ClassifierTrainer, \
    UnsupervisedTrainer, adversarial_wrapper, AdapterTrainer
from scvi.inference.annotation import compute_accuracy_rf, compute_accuracy_svc
from scvi.models import VAE, SCANVI, VAEC
from scvi.models.classifier import Classifier

use_cuda = True


def test_cortex():
    cortex_dataset = CortexDataset(save_path='tests/data/')
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, train_size=0.9, use_cuda=use_cuda)
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.train_set.ll()
    trainer_cortex_vae.train_set.differential_expression_stats()

    trainer_cortex_vae.corrupt_posteriors(corruption='binomial')
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.uncorrupt_posteriors()

    trainer_cortex_vae.train_set.imputation_benchmark(n_samples=1)

    svaec = SCANVI(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    trainer_cortex_svaec = JointSemiSupervisedTrainer(svaec, cortex_dataset,
                                                      n_labelled_samples_per_class=50,
                                                      use_cuda=use_cuda)
    trainer_cortex_svaec.train(n_epochs=1)
    trainer_cortex_svaec.labelled_set.accuracy()
    trainer_cortex_svaec.full_dataset.ll()

    svaec = SCANVI(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    trainer_cortex_svaec = AlternateSemiSupervisedTrainer(svaec, cortex_dataset,
                                                          n_labelled_samples_per_class=50,
                                                          use_cuda=use_cuda)
    trainer_cortex_svaec.train(n_epochs=1, lr=1e-2)
    trainer_cortex_svaec.unlabelled_set.accuracy()
    data_train, labels_train = trainer_cortex_svaec.labelled_set.raw_data()
    data_test, labels_test = trainer_cortex_svaec.unlabelled_set.raw_data()
    compute_accuracy_svc(data_train, labels_train, data_test, labels_test,
                         param_grid=[{'C': [1], 'kernel': ['linear']}])
    compute_accuracy_rf(data_train, labels_train, data_test, labels_test,
                        param_grid=[{'max_depth': [3], 'n_estimators': [10]}])

    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    cls_trainer = ClassifierTrainer(cls, cortex_dataset)
    cls_trainer.train(n_epochs=1)
    cls_trainer.train_set.accuracy()


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()
    synthetic_dataset.cell_types = np.array(['A', 'B', 'C'])
    svaec = SCANVI(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    trainer_synthetic_svaec = JointSemiSupervisedTrainer(svaec, synthetic_dataset, use_cuda=use_cuda)
    trainer_synthetic_svaec.train(n_epochs=1)
    trainer_synthetic_svaec.labelled_set.entropy_batch_mixing()
    trainer_synthetic_svaec.full_dataset.knn_purity(verbose=True)
    trainer_synthetic_svaec.labelled_set.show_t_sne(n_samples=50)
    trainer_synthetic_svaec.unlabelled_set.show_t_sne(n_samples=50, color_by='labels')
    trainer_synthetic_svaec.labelled_set.show_t_sne(n_samples=50, color_by='batches and labels')
    trainer_synthetic_svaec.labelled_set.clustering_scores()
    trainer_synthetic_svaec.labelled_set.clustering_scores(prediction_algorithm='gmm')
    trainer_synthetic_svaec.unlabelled_set.unsupervised_accuracy()
    trainer_synthetic_svaec.unlabelled_set.differential_expression_score('B', 'C', genes=['2', '4'])
    trainer_synthetic_svaec.unlabelled_set.differential_expression_table()


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    vaec = VAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    trainer_synthetic_vaec = JointSemiSupervisedTrainer(vaec, synthetic_dataset, use_cuda=use_cuda, frequency=1,
                                                        early_stopping_kwargs={'early_stopping_metric': 'll',
                                                                               'on': 'labelled_set',
                                                                               'save_best_state_metric': 'll'})
    trainer_synthetic_vaec = adversarial_wrapper(trainer_synthetic_vaec, warm_up=1)
    trainer_synthetic_vaec.train(n_epochs=2)


def test_fish_rna():
    gene_dataset_fish = SmfishDataset('tests/data/')
    gene_dataset_seq = CortexDataset(save_path='tests/data/',
                                     genes_fish=gene_dataset_fish.gene_names,
                                     genes_to_keep=[], additional_genes=50)
    benchamrk_fish_scrna(gene_dataset_seq, gene_dataset_fish)


def base_benchmark(gene_dataset):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=1)
    return trainer


def test_all_benchmarks():
    all_benchmarks(n_epochs=1)


def test_synthetic_3():
    gene_dataset = SyntheticDataset()
    trainer = base_benchmark(gene_dataset)
    adapter_trainer = AdapterTrainer(trainer.model, gene_dataset, trainer.train_set, frequency=1)
    adapter_trainer.train(n_path=1, n_epochs=1)


def test_brain_large():
    brain_large_dataset = BrainLargeDataset(subsample_size=128, save_path='tests/data/')
    base_benchmark(brain_large_dataset)


def test_retina():
    retina_dataset = RetinaDataset(save_path='tests/data/')
    base_benchmark(retina_dataset)


def test_cite_seq():
    pbmc_cite_seq_dataset = CiteSeqDataset(name='pbmc', save_path='tests/data/citeSeq/')
    base_benchmark(pbmc_cite_seq_dataset)


def test_brain_small():
    brain_small_dataset = BrainSmallDataset(save_path='tests/data/')
    base_benchmark(brain_small_dataset)


def test_hemato():
    hemato_dataset = HematoDataset(save_path='tests/data/HEMATO/')
    base_benchmark(hemato_dataset)


def test_loom():
    retina_dataset = LoomDataset("retina.loom", save_path='tests/data/')
    base_benchmark(retina_dataset)


def test_remote_loom():
    fish_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom",
                               save_path='tests/data/',
                               url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
    base_benchmark(fish_dataset)


def test_cortex_loom():
    cortex_dataset = LoomDataset("Cortex.loom",
                                 save_path='tests/data/')
    base_benchmark(cortex_dataset)


def test_anndata():
    ann_dataset = AnnDataset("TM_droplet_mat.h5ad", save_path='tests/data/')
    base_benchmark(ann_dataset)


def test_csv():
    csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='tests/data/', compression='gzip')
    base_benchmark(csv_dataset)


def test_cbmc():
    cbmc_dataset = CbmcDataset(save_path='tests/data/citeSeq/')
    trainer = base_benchmark(cbmc_dataset)
    trainer.train_set.nn_overlap_score(k=5)


def test_pbmc():
    pbmc_dataset = PbmcDataset(save_path='tests/data/')
    purified_pbmc_dataset = PurifiedPBMCDataset(save_path='tests/data/')  # all cells
    purified_t_cells = PurifiedPBMCDataset(save_path='tests/data/', filter_cell_types=range(6))  # only t-cells
    base_benchmark(pbmc_dataset)
    assert len(purified_t_cells.cell_types) == 6
    assert len(purified_pbmc_dataset.cell_types) == 10


def test_filter_and_concat_datasets():
    cortex_dataset_1 = CortexDataset(save_path='tests/data/')
    cortex_dataset_1.subsample_genes(subset_genes=np.arange(0, 10))
    cortex_dataset_1.filter_cell_types(["microglia", "oligodendrocytes"])
    cortex_dataset_2 = CortexDataset(save_path='tests/data/')
    cortex_dataset_2.subsample_genes(subset_genes=np.arange(10, 20))
    cortex_dataset_2.filter_cell_types(["endothelial-mural", "interneurons", "microglia", "oligodendrocytes"])
    cortex_dataset_2.filter_cell_types([2, 0])
    cortex_dataset_merged = GeneExpressionDataset.concat_datasets(cortex_dataset_1, cortex_dataset_2)
    assert cortex_dataset_merged.nb_genes == 20

    synthetic_dataset_1 = SyntheticDataset(n_batches=2, n_labels=5)
    synthetic_dataset_2 = SyntheticDataset(n_batches=3, n_labels=3)
    synthetic_merged_1 = GeneExpressionDataset.concat_datasets(synthetic_dataset_1, synthetic_dataset_2)
    assert synthetic_merged_1.n_batches == 5
    assert synthetic_merged_1.n_labels == 5

    synthetic_merged_2 = GeneExpressionDataset.concat_datasets(synthetic_dataset_1, synthetic_dataset_2,
                                                               shared_labels=False)
    assert synthetic_merged_2.n_batches == 5
    assert synthetic_merged_2.n_labels == 8

    synthetic_dataset_1.filter_cell_types([0, 1, 2, 3])
    assert synthetic_dataset_1.n_labels == 4

    synthetic_dataset_1.subsample_cells(50)
    assert len(synthetic_dataset_1) == 50

    synthetic_dataset_3 = SyntheticDataset(n_labels=6)
    synthetic_dataset_3.cell_types = np.arange(6).astype(np.str)
    synthetic_dataset_3.map_cell_types({"2": "9", ("4", "3"): "8"})


def test_seqfish():
    seqfish_dataset = SeqfishDataset(save_path='tests/data/')
    base_benchmark(seqfish_dataset)


def test_breast_cancer():
    breast_cancer_dataset = BreastCancerDataset(save_path='tests/data/')
    base_benchmark(breast_cancer_dataset)


def test_mouseob():
    mouseob_dataset = MouseOBDataset(save_path='tests/data/')
    base_benchmark(mouseob_dataset)


def test_particular_benchmark():
    synthetic_dataset = SyntheticDataset()
    benchmark(synthetic_dataset, n_epochs=1, use_cuda=False)


def test_nb_not_zinb():
    synthetic_dataset = SyntheticDataset()
    svaec = SCANVI(synthetic_dataset.nb_genes,
                   synthetic_dataset.n_batches,
                   synthetic_dataset.n_labels,
                   labels_groups=[0, 0, 1],
                   reconstruction_loss="nb")
    trainer_synthetic_svaec = JointSemiSupervisedTrainer(svaec, synthetic_dataset, use_cuda=use_cuda)
    trainer_synthetic_svaec.train(n_epochs=1)
