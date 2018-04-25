#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

from run_benchmarks import run_benchmarks
from scvi.dataset import load_datasets


def test_benchmark():
    gene_dataset = load_datasets("synthetic")
    run_benchmarks(gene_dataset, n_epochs=1, use_batches=True)


def test_cortex():
    gene_dataset = load_datasets("cortex")
    run_benchmarks(gene_dataset, n_epochs=1, semi_supervised=True)


def test_brain_large():
    gene_dataset = load_datasets("brain_large", unit_test=True)
    run_benchmarks(gene_dataset, n_epochs=1)
