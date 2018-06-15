#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""
from scvi.benchmark import run_benchmarks
from scvi.benchmark import run_benchmarks_classification
from scvi.models import VAEC, VAE, SVAEC
from scvi.dataset import load_datasets


def test_synthetic_1():
    synthetic_dataset = load_datasets("synthetic")
    run_benchmarks(synthetic_dataset, n_epochs=1, use_batches=True, model=VAE)
    run_benchmarks_classification(synthetic_dataset, n_epochs=1, n_epochs_classifier=1)


def test_synthetic_2():
    synthetic_dataset = load_datasets("synthetic")
    run_benchmarks(synthetic_dataset, n_epochs=1, model=SVAEC, benchmark=True)


def test_cortex():
    cortex_dataset = load_datasets("cortex")
    run_benchmarks(cortex_dataset, n_epochs=1, model=VAEC)


def test_brain_large():
    brain_large_dataset = load_datasets("brain_large", save_path='tests/data/')
    run_benchmarks(brain_large_dataset, n_epochs=1, use_batches=False, tt_split=0.5)


def test_retina():
    retina_dataset = load_datasets("retina", save_path='tests/data/')
    run_benchmarks(retina_dataset, n_epochs=1, show_batch_mixing=False)


def test_cbmc():
    cbmc_dataset = load_datasets("cbmc", save_path='tests/data/')
    run_benchmarks(cbmc_dataset, n_epochs=1, show_batch_mixing=False)


def test_brain_small():
    brain_small_dataset = load_datasets("brain_small", save_path='tests/data/')
    run_benchmarks(brain_small_dataset, n_epochs=1, show_batch_mixing=False)


def test_hemato():
    hemato_dataset = load_datasets("hemato", save_path='tests/data/HEMATO/')
    run_benchmarks(hemato_dataset, n_epochs=1, show_batch_mixing=False)


def test_pbmc():
    pbmc_dataset = load_datasets("pbmc", save_path='tests/data/PBMC/')
    run_benchmarks(pbmc_dataset, n_epochs=1, show_batch_mixing=False)


def test_loom():
    retina_dataset = load_datasets("retina.loom", save_path='tests/data/',
                                   batch_col='Batch_id', cluster_col='Clusters')
    run_benchmarks(retina_dataset, n_epochs=1, show_batch_mixing=False)


def test_remote_loom():
    fish_dataset = load_datasets("osmFISH_SScortex_mouse_all_cell.loom", gene_row='Gene',
                                 cluster_col='ClusterID', save_path='data/',
                                 url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
    run_benchmarks(fish_dataset, n_epochs=10, show_batch_mixing=False)
