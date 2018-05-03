#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""
from scvi.benchmark import run_benchmarks_classification
from run_benchmarks import run_benchmarks
from scvi.dataset import load_datasets
from scvi.models import VAEC, VAE, SVAEC


def test_synthetic_1():
    gene_dataset = load_datasets("synthetic")
    run_benchmarks(gene_dataset, n_epochs=1, use_batches=True, model=VAE)
    run_benchmarks_classification(gene_dataset, n_epochs=1)


def test_synthetic_2():
    gene_dataset = load_datasets("synthetic")
    run_benchmarks(gene_dataset, n_epochs=1, model=SVAEC, benchmark=True)


def test_cortex():
    gene_dataset = load_datasets("cortex")
    run_benchmarks(gene_dataset, n_epochs=1, model=VAEC)


def test_brain_large():
    gene_dataset = load_datasets("brain_large", unit_test=True)
    run_benchmarks(gene_dataset, n_epochs=1)


def test_retina():
    gene_dataset = load_datasets("retina", unit_test=True)
    run_benchmarks(gene_dataset, n_epochs=1, show_batch_mixing=False)
