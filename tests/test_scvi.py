#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""
import numpy as np

from scvi.benchmark import run_benchmarks, run_benchmarks_classification
from scvi.dataset import BrainLargeDataset, CortexDataset, SyntheticDataset, \
    RetinaDataset, BrainSmallDataset, HematoDataset, LoomDataset, AnnDataset, CsvDataset, \
    CiteSeqDataset, GeneExpressionDataset
from scvi.models import VAEC, VAE, SVAEC


def test_synthetic_1():
    synthetic_dataset = SyntheticDataset()
    run_benchmarks(synthetic_dataset, n_epochs=1, use_batches=True, model=VAE)
    run_benchmarks_classification(synthetic_dataset, n_epochs=1, n_epochs_classifier=1)


def test_synthetic_2():
    synthetic_dataset = SyntheticDataset()
    run_benchmarks(synthetic_dataset, n_epochs=1, model=SVAEC, benchmark=True)


def test_cortex():
    cortex_dataset = CortexDataset()
    run_benchmarks(cortex_dataset, n_epochs=1, model=VAEC)


def test_brain_large():
    brain_large_dataset = BrainLargeDataset(subsample_size=128, save_path='tests/data/')
    run_benchmarks(brain_large_dataset, n_epochs=1, use_batches=False, tt_split=0.5)


def test_retina():
    retina_dataset = RetinaDataset(save_path='tests/data/')
    run_benchmarks(retina_dataset, n_epochs=1, show_batch_mixing=False)


def test_cite_seq():
    pbmc_cite_seq_dataset = CiteSeqDataset(name='pbmc', save_path='tests/data/citeSeq/')
    run_benchmarks(pbmc_cite_seq_dataset, n_epochs=1, show_batch_mixing=False)


def test_brain_small():
    brain_small_dataset = BrainSmallDataset(save_path='tests/data/')
    run_benchmarks(brain_small_dataset, n_epochs=1, show_batch_mixing=False)


def test_hemato():
    hemato_dataset = HematoDataset(save_path='tests/data/HEMATO/')
    run_benchmarks(hemato_dataset, n_epochs=1, show_batch_mixing=False)


def test_loom():
    retina_dataset = LoomDataset("retina.loom", save_path='tests/data/')
    run_benchmarks(retina_dataset, n_epochs=1, show_batch_mixing=False)


def test_remote_loom():
    fish_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom",
                               save_path='data/',
                               url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
    run_benchmarks(fish_dataset, n_epochs=10, show_batch_mixing=False)


def test_cortex_loom():
    cortex_dataset = LoomDataset("Cortex.loom",
                                 save_path='tests/data/')
    run_benchmarks(cortex_dataset, n_epochs=1, show_batch_mixing=False)


def test_anndata():
    ann_dataset = AnnDataset("test.h5ad", save_path='tests/data/')
    run_benchmarks(ann_dataset, n_epochs=1, show_batch_mixing=False)


def test_csv():
    csv_dataet = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='tests/data/', compression='gzip')
    run_benchmarks(csv_dataet, n_epochs=1, show_batch_mixing=False)


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
