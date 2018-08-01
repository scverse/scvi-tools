#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np

from scvi.benchmark import all_benchmarks, benchmark
from scvi.dataset import BrainLargeDataset, CortexDataset, RetinaDataset, BrainSmallDataset, HematoDataset, \
    LoomDataset, AnnDataset, CsvDataset, CiteSeqDataset, CbmcDataset, PbmcDataset, SyntheticDataset, \
    SeqfishDataset, SmfishDataset, BreastCancerDataset, MouseOBDataset, \
    GeneExpressionDataset, PurifiedPBMCDataset
from smFISH.dataset import CortexDatasetCustom, SmfishDatasetCustom
from smFISH.metrics.visualisation import show_spatial_expression, show_gene_exp, \
    show_mixing, show_cell_types, compare_cell_types
from scvi.inference import JointSemiSupervisedVariationalInference, AlternateSemiSupervisedVariationalInference, \
    ClassifierInference, VariationalInference, adversarial_wrapper, mmd_wrapper, VariationalInferenceFish
from scvi.metrics.adapt_encoder import adapt_encoder
from scvi.models import VAE, SVAEC, VAEC, VAEF
from scvi.models.classifier import Classifier
from scvi.metrics.clustering import get_data, get_common_t_sne

use_cuda = True


def test_cortex():
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.1, use_cuda=use_cuda)
    infer_cortex_vae.train(n_epochs=1)
    infer_cortex_vae.ll('train')
    infer_cortex_vae.differential_expression_stats('train')
    infer_cortex_vae.differential_expression('test')
    infer_cortex_vae.imputation('train', corruption='uniform')
    infer_cortex_vae.imputation('test', n_samples=2, corruption='binomial')

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    infer_cortex_svaec = JointSemiSupervisedVariationalInference(svaec, cortex_dataset,
                                                                 n_labelled_samples_per_class=50,
                                                                 use_cuda=use_cuda)
    infer_cortex_svaec.train(n_epochs=1)
    infer_cortex_svaec.accuracy('labelled')
    infer_cortex_svaec.ll('all')

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels,
                  logreg_classifier=True)
    infer_cortex_svaec = AlternateSemiSupervisedVariationalInference(svaec, cortex_dataset,
                                                                     n_labelled_samples_per_class=50,
                                                                     use_cuda=use_cuda)
    infer_cortex_svaec.train(n_epochs=1, lr=1e-2)
    infer_cortex_svaec.accuracy('unlabelled')
    infer_cortex_svaec.svc_rf(unit_test=True)

    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    infer_cls = ClassifierInference(cls, cortex_dataset)
    infer_cls.train(n_epochs=1)
    infer_cls.accuracy('train')


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()
    svaec = SVAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    infer_synthetic_svaec = JointSemiSupervisedVariationalInference(svaec, synthetic_dataset, use_cuda=use_cuda)
    infer_synthetic_svaec.train(n_epochs=1)
    infer_synthetic_svaec.entropy_batch_mixing('labelled')
    infer_synthetic_svaec.show_t_sne('labelled', n_samples=50)
    infer_synthetic_svaec.show_t_sne('unlabelled', n_samples=50, color_by='labels')
    infer_synthetic_svaec.show_t_sne('labelled', n_samples=50, color_by='batches and labels')
    infer_synthetic_svaec.clustering_scores('labelled')
    infer_synthetic_svaec.clustering_scores('labelled', prediction_algorithm='gmm')
    infer_synthetic_svaec.unsupervised_accuracy('unlabelled')


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    vaec = VAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    infer_synthetic_vaec = JointSemiSupervisedVariationalInference(vaec, synthetic_dataset, use_cuda=use_cuda,
                                                                   early_stopping_metric='ll', frequency=1,
                                                                   save_best_state_metric='accuracy', on='labelled')
    infer_synthetic_vaec = adversarial_wrapper(infer_synthetic_vaec, warm_up=5)
    infer_synthetic_vaec = mmd_wrapper(infer_synthetic_vaec, warm_up=15)
    infer_synthetic_vaec.train(n_epochs=20)
    infer_synthetic_vaec.svc_rf(unit_test=True)


def test_fish_rna():
    gene_dataset_fish = SmfishDatasetCustom()
    gene_names = gene_dataset_fish.gene_names
    indexes_to_keep = np.arange(len(gene_names))
    gene_dataset_seq = CortexDatasetCustom(genes_fish=gene_dataset_fish.gene_names,
                                           genes_to_keep=[], additional_genes=50)
    vae = VAEF(gene_dataset_seq.nb_genes, indexes_to_keep, n_layers_decoder=2, n_latent=6,
               n_layers=2, n_hidden=256, reconstruction_loss='nb', dropout_rate=0.3, n_labels=7, n_batch=0,
               model_library=False)
    infer = VariationalInferenceFish(vae, gene_dataset_seq, gene_dataset_fish, train_size=0.9, verbose=True,
                                     frequency=5, weight_decay=0.35, n_epochs_even=100, n_epochs_kl=1000,
                                     cl_ratio=0, n_epochs_cl=100)
    infer.train(n_epochs=1, lr=0.0008)
    data_loader_fish = infer.data_loaders['train_fish']
    data_loader_seq = infer.data_loaders['train_seq']
    latent_seq, _, labels_seq, expected_frequencies_seq, values_seq = get_data(vae, data_loader_seq, mode="scRNA")
    latent_fish, _, labels_fish, expected_frequencies_fish, values_fish, \
        x_coords, y_coords = get_data(vae, data_loader_fish, mode="smFISH")
    t_sne_seq, t_sne_fish, idx_t_sne_seq, idx_t_sne_fish = get_common_t_sne(latent_seq, latent_fish, n_samples=10)
    show_cell_types(t_sne_seq, labels_seq[idx_t_sne_seq], t_sne_fish, labels_fish[idx_t_sne_fish])
    show_mixing(t_sne_seq, t_sne_fish)
    show_gene_exp(t_sne_seq, expected_frequencies_seq[idx_t_sne_seq, 0], labels=labels_seq[idx_t_sne_seq],
                  title="", gene_name="")
    compare_cell_types(t_sne_fish, labels_fish[idx_t_sne_fish], labels_fish[idx_t_sne_fish])
    show_spatial_expression(x_coords, y_coords,
                            expected_frequencies_fish[:, 0], labels=labels_fish,
                            title="", gene_name="gad2")


def base_benchmark(gene_dataset):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    infer = VariationalInference(vae, gene_dataset, train_size=0.5, use_cuda=use_cuda)
    infer.train(n_epochs=1)
    return infer


def test_all_benchmarks():
    all_benchmarks(n_epochs=1, unit_test=True)


def test_synthetic_3():
    infer = base_benchmark(SyntheticDataset())
    adapt_encoder(infer, n_path=1, n_epochs=1, frequency=1)


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
    ann_dataset = AnnDataset("test.h5ad", save_path='tests/data/')
    base_benchmark(ann_dataset)


def test_csv():
    csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='tests/data/', compression='gzip')
    base_benchmark(csv_dataset)


def test_cbmc():
    cbmc_dataset = CbmcDataset(save_path='tests/data/citeSeq/')
    infer = base_benchmark(cbmc_dataset)
    infer.nn_overlap_score(k=5)


def test_pbmc():
    pbmc_dataset = PbmcDataset(save_path='tests/data/')
    purified_pbmc_dataset = PurifiedPBMCDataset(save_path='tests/data/')  # all cells
    purified_t_cells = PurifiedPBMCDataset(save_path='tests/data/', filter_cell_types=range(6))  # only t-cells
    base_benchmark(pbmc_dataset)
    assert len(purified_t_cells.cell_types) == 6
    assert len(purified_pbmc_dataset.cell_types) == 10


def test_filter_and_concat_datasets():
    cortex_dataset_1 = CortexDataset()
    cortex_dataset_1.subsample_genes(subset_genes=np.arange(0, 300))
    cortex_dataset_1.filter_cell_types(["microglia", "oligodendrocytes"])
    cortex_dataset_2 = CortexDataset()
    cortex_dataset_2.subsample_genes(subset_genes=np.arange(100, 400))
    cortex_dataset_2.filter_cell_types(["endothelial-mural", "interneurons", "microglia", "oligodendrocytes"])
    cortex_dataset_2.filter_cell_types([2, 0])
    cortex_dataset_merged = GeneExpressionDataset.concat_datasets(cortex_dataset_1, cortex_dataset_2)
    assert cortex_dataset_merged.nb_genes == 200

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


def test_seqfish():
    seqfish_dataset = SeqfishDataset(save_path='tests/data/')
    base_benchmark(seqfish_dataset)


def test_breast_cancer():
    breast_cancer_dataset = BreastCancerDataset(save_path='tests/data/')
    base_benchmark(breast_cancer_dataset)


def test_mouseob():
    mouseob_dataset = MouseOBDataset(save_path='tests/data/')
    base_benchmark(mouseob_dataset)


def test_smfish():
    smfish_dataset = SmfishDataset(save_path='tests/data/')
    base_benchmark(smfish_dataset)


def test_particular_benchmark():
    synthetic_dataset = SyntheticDataset()
    benchmark(synthetic_dataset, n_epochs=1, use_cuda=False)


def test_nb_not_zinb():
    synthetic_dataset = SyntheticDataset()
    svaec = SVAEC(synthetic_dataset.nb_genes,
                  synthetic_dataset.n_batches,
                  synthetic_dataset.n_labels,
                  labels_groups=[0, 0, 1],
                  reconstruction_loss="nb")
    infer_synthetic_svaec = JointSemiSupervisedVariationalInference(svaec, synthetic_dataset, use_cuda=use_cuda)
    infer_synthetic_svaec.train(n_epochs=1)
