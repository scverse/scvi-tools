#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np

from scvi.benchmark import all_benchmarks, benchmark, benchmark_fish_scrna
from scvi.dataset import BrainLargeDataset, CortexDataset, RetinaDataset, BrainSmallDataset, HematoDataset, \
    LoomDataset, AnnDataset, CsvDataset, CiteSeqDataset, CbmcDataset, PbmcDataset, SyntheticDataset, \
    SeqfishDataset, SmfishDataset, BreastCancerDataset, MouseOBDataset, \
    GeneExpressionDataset, PurifiedPBMCDataset
from scvi.inference import JointSemiSupervisedTrainer, AlternateSemiSupervisedTrainer, ClassifierTrainer, \
    UnsupervisedTrainer, AdapterTrainer
from scvi.inference.annotation import compute_accuracy_rf, compute_accuracy_svc
from scvi.models import VAE, SCANVI, VAEC
from scvi.models.classifier import Classifier
import anndata
import os.path

use_cuda = True


def test_cortex(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    trainer_cortex_vae = UnsupervisedTrainer(vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda)
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.train_set.ll()
    trainer_cortex_vae.train_set.differential_expression_stats()

    trainer_cortex_vae.corrupt_posteriors(corruption='binomial')
    trainer_cortex_vae.corrupt_posteriors()
    trainer_cortex_vae.train(n_epochs=1)
    trainer_cortex_vae.uncorrupt_posteriors()

    trainer_cortex_vae.train_set.imputation_benchmark(n_samples=1, show_plot=False,
                                                      title_plot='imputation', save_path=save_path)

    svaec = SCANVI(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    trainer_cortex_svaec = JointSemiSupervisedTrainer(svaec, cortex_dataset,
                                                      n_labelled_samples_per_class=3,
                                                      use_cuda=use_cuda)
    trainer_cortex_svaec.train(n_epochs=1)
    trainer_cortex_svaec.labelled_set.accuracy()
    trainer_cortex_svaec.full_dataset.ll()

    svaec = SCANVI(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    trainer_cortex_svaec = AlternateSemiSupervisedTrainer(svaec, cortex_dataset,
                                                          n_labelled_samples_per_class=3,
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
    trainer_synthetic_svaec.labelled_set.show_t_sne(n_samples=5)
    trainer_synthetic_svaec.unlabelled_set.show_t_sne(n_samples=5, color_by='labels')
    trainer_synthetic_svaec.labelled_set.show_t_sne(n_samples=5, color_by='batches and labels')
    trainer_synthetic_svaec.labelled_set.clustering_scores()
    trainer_synthetic_svaec.labelled_set.clustering_scores(prediction_algorithm='gmm')
    trainer_synthetic_svaec.unlabelled_set.unsupervised_classification_accuracy()
    trainer_synthetic_svaec.unlabelled_set.differential_expression_score(synthetic_dataset.labels.ravel() == 1,
                                                                         synthetic_dataset.labels.ravel() == 2,
                                                                         genes=['2', '4'], n_samples=2,
                                                                         M_permutation=10)
    trainer_synthetic_svaec.unlabelled_set.one_vs_all_degenes(n_samples=2, M_permutation=10)


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    vaec = VAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    trainer_synthetic_vaec = JointSemiSupervisedTrainer(vaec, synthetic_dataset, use_cuda=use_cuda, frequency=1,
                                                        early_stopping_kwargs={'early_stopping_metric': 'll',
                                                                               'on': 'labelled_set',
                                                                               'save_best_state_metric': 'll'})
    trainer_synthetic_vaec.train(n_epochs=2)


def test_fish_rna(save_path):
    gene_dataset_fish = SmfishDataset(save_path)
    gene_dataset_seq = CortexDataset(save_path=save_path,
                                     genes_fish=gene_dataset_fish.gene_names,
                                     genes_to_keep=[], additional_genes=50)
    benchmark_fish_scrna(gene_dataset_seq, gene_dataset_fish)


def base_benchmark(gene_dataset):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=1)
    return trainer


def test_all_benchmarks(save_path):
    all_benchmarks(n_epochs=1, save_path=save_path, show_plot=False)


def test_synthetic_3():
    gene_dataset = SyntheticDataset()
    trainer = base_benchmark(gene_dataset)
    adapter_trainer = AdapterTrainer(trainer.model, gene_dataset, trainer.train_set, frequency=1)
    adapter_trainer.train(n_path=1, n_epochs=1)


def test_brain_large(save_path):
    brain_large_dataset = BrainLargeDataset(subsample_size=128, save_path=save_path)
    base_benchmark(brain_large_dataset)


def test_retina(save_path):
    retina_dataset = RetinaDataset(save_path=save_path)
    base_benchmark(retina_dataset)


def test_cite_seq(save_path):
    pbmc_cite_seq_dataset = CiteSeqDataset(name='pbmc', save_path=os.path.join(save_path, 'citeSeq/'))
    base_benchmark(pbmc_cite_seq_dataset)


def test_brain_small(save_path):
    brain_small_dataset = BrainSmallDataset(save_path=save_path)
    base_benchmark(brain_small_dataset)


def test_hemato(save_path):
    hemato_dataset = HematoDataset(save_path=os.path.join(save_path, 'HEMATO/'))
    base_benchmark(hemato_dataset)


def test_loom(save_path):
    retina_dataset = LoomDataset("retina.loom", save_path=save_path)
    base_benchmark(retina_dataset)


def test_remote_loom(save_path):
    fish_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom",
                               save_path=save_path,
                               url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
    base_benchmark(fish_dataset)


def test_cortex_loom(save_path):
    cortex_dataset = LoomDataset("Cortex.loom",
                                 save_path=save_path)
    base_benchmark(cortex_dataset)


def test_anndata(save_path):
    ann_dataset = AnnDataset("TM_droplet_mat.h5ad", save_path=save_path)
    base_benchmark(ann_dataset)
    AnnDataset(anndata.AnnData(np.random.randint(1, 10, (10, 10))))


def test_csv(save_path):
    csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path=save_path, compression='gzip')
    base_benchmark(csv_dataset)


def test_cbmc(save_path):
    cbmc_dataset = CbmcDataset(save_path=os.path.join(save_path, 'citeSeq/'))
    trainer = base_benchmark(cbmc_dataset)
    trainer.train_set.nn_overlap_score(k=5)


def test_pbmc(save_path):
    pbmc_dataset = PbmcDataset(save_path=save_path)
    purified_pbmc_dataset = PurifiedPBMCDataset(save_path=save_path)  # all cells
    purified_t_cells = PurifiedPBMCDataset(save_path=save_path, filter_cell_types=range(6))  # only t-cells
    base_benchmark(pbmc_dataset)
    assert len(purified_t_cells.cell_types) == 6
    assert len(purified_pbmc_dataset.cell_types) == 10


def test_filter_and_concat_datasets(save_path):
    cortex_dataset_1 = CortexDataset(save_path=save_path)
    cortex_dataset_1.subsample_genes(subset_genes=np.arange(0, 3))
    cortex_dataset_1.filter_cell_types(["microglia", "oligodendrocytes"])
    cortex_dataset_2 = CortexDataset(save_path=save_path)
    cortex_dataset_2.subsample_genes(subset_genes=np.arange(1, 4))
    cortex_dataset_2.filter_cell_types(["endothelial-mural", "interneurons", "microglia", "oligodendrocytes"])
    cortex_dataset_2.filter_cell_types([2, 0])
    cortex_dataset_merged = GeneExpressionDataset.concat_datasets(cortex_dataset_1, cortex_dataset_2)
    assert cortex_dataset_merged.nb_genes == 2

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


def test_seqfish(save_path):
    seqfish_dataset = SeqfishDataset(save_path=save_path)
    base_benchmark(seqfish_dataset)


def test_breast_cancer(save_path):
    breast_cancer_dataset = BreastCancerDataset(save_path=save_path)
    base_benchmark(breast_cancer_dataset)


def test_mouseob(save_path):
    mouseob_dataset = MouseOBDataset(save_path=save_path)
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


def test_classifier_accuracy(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    cls_trainer = ClassifierTrainer(cls, cortex_dataset, metrics_to_monitor=['accuracy'], frequency=1,
                                    early_stopping_kwargs={'early_stopping_metric': 'accuracy',
                                                           'save_best_state_metric': 'accuracy'})
    cls_trainer.train(n_epochs=2)
    cls_trainer.train_set.accuracy()


def test_sampling_zl(save_path):
    cortex_dataset = CortexDataset(save_path=save_path)
    cortex_vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    trainer_cortex_vae = UnsupervisedTrainer(cortex_vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda)
    trainer_cortex_vae.train(n_epochs=2)

    cortex_cls = Classifier((cortex_vae.n_latent + 1), n_labels=cortex_dataset.n_labels)
    trainer_cortex_cls = ClassifierTrainer(cortex_cls, cortex_dataset,
                                           sampling_model=cortex_vae, sampling_zl=True)
    trainer_cortex_cls.train(n_epochs=2)
    trainer_cortex_cls.test_set.accuracy()
