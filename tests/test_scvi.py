#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np

from scvi.benchmark import UnsupervisedInference, SemiSupervisedInference
from scvi.dataset import BrainLargeDataset, CortexDataset, RetinaDataset, BrainSmallDataset, HematoDataset, \
    LoomDataset, AnnDataset, CsvDataset, CiteSeqDataset, CbmcDataset, PbmcDataset, SyntheticDataset, \
    GeneExpressionDataset
from scvi.inference import JointSemiSupervisedVariationalInference, AlternateSemiSupervisedVariationalInference, \
    ClassifierInference, VariationalInference
from scvi.models import VAE, VAEC, SVAEC
from scvi.models.classifier import Classifier

use_cuda = True


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()

    vae = VAE(synthetic_dataset.nb_genes, n_batch=synthetic_dataset.n_batches)
    infer = UnsupervisedInference(synthetic_dataset, vae)
    infer.train(n_epochs=1, benchmark=False, verbose=True)
    infer.all(unit_test=True)

    svaec = SVAEC(synthetic_dataset.nb_genes, n_batch=synthetic_dataset.n_batches, n_labels=synthetic_dataset.n_labels,
                  logreg_classifier=True)
    infer = SemiSupervisedInference(synthetic_dataset, svaec)
    infer.train(n_epochs=1, benchmark=False, mode="alternately", verbose=True)
    infer.all(unit_test=True)

    svaec = SVAEC(synthetic_dataset.nb_genes, n_batch=synthetic_dataset.n_batches, n_labels=synthetic_dataset.n_labels,
                  logreg_classifier=False)
    infer = SemiSupervisedInference(synthetic_dataset, svaec)
    infer.train(n_epochs=1, benchmark=False, mode="jointly", verbose=True)
    infer.all(unit_test=True)

    vaec = VAEC(synthetic_dataset.nb_genes, n_batch=synthetic_dataset.n_batches, n_labels=synthetic_dataset.n_labels)
    infer = SemiSupervisedInference(synthetic_dataset, vaec)
    infer.train(n_epochs=1, benchmark=False, mode="jointly", verbose=True)
    infer.all(unit_test=True)

def test_cortex():
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_labels, n_batch=cortex_dataset.n_batches)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.1)
    infer_cortex_vae.fit(n_epochs=1)
    infer_cortex_vae.ll('train')
    infer_cortex_vae.de_stats('train')
    infer_cortex_vae.de_stats('test')
    infer_cortex_vae.imputation_stats('test', rate=0.5)

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_labels, n_batch=cortex_dataset.n_batches)
    infer_cortex_svaec = JointSemiSupervisedVariationalInference(svaec, cortex_dataset, n_labelled_samples_per_class=50)
    infer_cortex_svaec.fit(n_epochs=1)
    infer_cortex_svaec.accuracy('labelled')
    infer_cortex_svaec.ll('all')

    svaec = SVAEC(cortex_dataset.nb_genes, cortex_dataset.n_labels, n_batch=cortex_dataset.n_batches,
                  logreg_classifier=True)
    infer_cortex_svaec = AlternateSemiSupervisedVariationalInference(svaec, cortex_dataset,
                                                                     n_labelled_samples_per_class=50)
    infer_cortex_svaec.fit(n_epochs=1, lr=1e-2)
    infer_cortex_svaec.accuracy('unlabelled')

    cls = Classifier(cortex_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
    infer_cls = ClassifierInference(cls, cortex_dataset)
    infer_cls.fit(n_epochs=1)
    infer_cls.accuracy('train')


def test_brain_large():
    brain_large_dataset = BrainLargeDataset(subsample_size=128, save_path='tests/data/')
    run_benchmarks(brain_large_dataset, use_batches=False, train_size=0.5, use_cuda=use_cuda)


def test_retina():
    retina_dataset = RetinaDataset(save_path='tests/data/')
    run_benchmarks(retina_dataset, use_cuda=use_cuda)


def test_cite_seq():
    pbmc_cite_seq_dataset = CiteSeqDataset(name='pbmc', save_path='tests/data/citeSeq/')
    run_benchmarks(pbmc_cite_seq_dataset, use_cuda=use_cuda)


def test_brain_small():
    brain_small_dataset = BrainSmallDataset(save_path='tests/data/')
    run_benchmarks(brain_small_dataset, use_cuda=use_cuda)


def test_hemato():
    hemato_dataset = HematoDataset(save_path='tests/data/HEMATO/')
    run_benchmarks(hemato_dataset, use_cuda=use_cuda)


def test_loom():
    retina_dataset = LoomDataset("retina.loom", save_path='tests/data/')
    run_benchmarks(retina_dataset, use_cuda=use_cuda)


def test_remote_loom():
    fish_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom",
                               save_path='tests/data/',
                               url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
    run_benchmarks(fish_dataset, use_cuda=use_cuda)


def test_cortex_loom():
    cortex_dataset = LoomDataset("Cortex.loom",
                                 save_path='tests/data/')
    run_benchmarks(cortex_dataset, use_cuda=use_cuda)


def test_anndata():
    ann_dataset = AnnDataset("test.h5ad", save_path='tests/data/')
    run_benchmarks(ann_dataset, use_cuda=use_cuda)


def test_csv():
    csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='tests/data/', compression='gzip')
    run_benchmarks(csv_dataset, use_cuda=use_cuda)


def test_cbmc():
    cbmc_dataset = CbmcDataset(save_path='tests/data/citeSeq/')
    run_benchmarks(cbmc_dataset, use_cuda=use_cuda)


def test_pbmc():
    pbmc_dataset = PbmcDataset(save_path='tests/data/')
    run_benchmarks(pbmc_dataset, use_cuda=use_cuda)


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
