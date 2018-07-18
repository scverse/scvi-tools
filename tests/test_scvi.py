#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np
import torch

from scvi.benchmark import all_benchmarks
from scvi.dataset import BrainLargeDataset, CortexDataset, RetinaDataset, BrainSmallDataset, HematoDataset, \
    LoomDataset, AnnDataset, CsvDataset, CiteSeqDataset, CbmcDataset, PbmcDataset, SyntheticDataset, \
    SeqfishDataset, SmfishDataset, BreastCancerDataset, MouseOBDataset, \
    GeneExpressionDataset
from scvi.inference import JointSemiSupervisedVariationalInference, AlternateSemiSupervisedVariationalInference, \
    ClassifierInference, VariationalInference
from scvi.metrics.adapt_encoder import adapt_encoder
from scvi.models import VAE, SVAEC, VAEC
from scvi.models.classifier import Classifier

use_cuda = True


def test_cortex():
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_batches)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.1, use_cuda=use_cuda)
    infer_cortex_vae.fit(n_epochs=1)
    infer_cortex_vae.ll('train')
    infer_cortex_vae.differential_expression_stats('train')
    infer_cortex_vae.differential_expression('test')
    infer_cortex_vae.imputation_errors('test', rate=0.5)

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels)
    infer_cortex_svaec = JointSemiSupervisedVariationalInference(svaec, cortex_dataset,
                                                                 n_labelled_samples_per_class=50,
                                                                 use_cuda=use_cuda)
    infer_cortex_svaec.fit(n_epochs=1)
    infer_cortex_svaec.accuracy('labelled')
    infer_cortex_svaec.ll('all')

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_batches, cortex_dataset.n_labels,
                  logreg_classifier=True)
    infer_cortex_svaec = AlternateSemiSupervisedVariationalInference(svaec, cortex_dataset,
                                                                     n_labelled_samples_per_class=50,
                                                                     use_cuda=use_cuda)
    infer_cortex_svaec.fit(n_epochs=1, lr=1e-2)
    infer_cortex_svaec.accuracy('unlabelled')
    infer_cortex_svaec.svc_rf(unit_test=True)

    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    infer_cls = ClassifierInference(cls, cortex_dataset)
    infer_cls.fit(n_epochs=1)
    infer_cls.accuracy('train')


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()
    svaec = SVAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    infer_synthetic_svaec = JointSemiSupervisedVariationalInference(svaec, synthetic_dataset, use_cuda=use_cuda)
    infer_synthetic_svaec.fit(n_epochs=1)
    with torch.no_grad():
        infer_synthetic_svaec.entropy_batch_mixing('labelled')
        infer_synthetic_svaec.show_t_sne('labelled', n_samples=50)
        infer_synthetic_svaec.show_t_sne('unlabelled', n_samples=50, color_by='labels')
        infer_synthetic_svaec.show_t_sne('labelled', n_samples=50, color_by='batches and labels')


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    vaec = VAEC(synthetic_dataset.nb_genes, synthetic_dataset.n_batches, synthetic_dataset.n_labels)
    infer_synthetic_vaec = JointSemiSupervisedVariationalInference(vaec, synthetic_dataset, use_cuda=use_cuda,
                                                                   early_stopping_metric='ll', frequency=1,
                                                                   save_best_state_metric='accuracy', on='labelled')
    infer_synthetic_vaec.fit(n_epochs=20)
    infer_synthetic_vaec.svc_rf(unit_test=True)


def base_benchmark(gene_dataset):
    vae = VAE(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels)
    infer = VariationalInference(vae, gene_dataset, train_size=0.5, use_cuda=use_cuda)
    infer.fit(n_epochs=1)
    return infer


def test_all_benchmarks():
    with torch.no_grad():
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
    base_benchmark(cbmc_dataset)


def test_pbmc():
    pbmc_dataset = PbmcDataset(save_path='tests/data/')
    base_benchmark(pbmc_dataset)


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
