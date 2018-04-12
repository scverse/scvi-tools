#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

from run_benchmarks import run_benchmarks
from scvi.dataset import load_datasets


def test_benchmark():
    gene_dataset_train, gene_dataset_test = load_datasets("synthetic")
    run_benchmarks(gene_dataset_train, gene_dataset_test, n_epochs=1)


def test_cortex():
    gene_dataset_train, gene_dataset_test = load_datasets("cortex")
    run_benchmarks(gene_dataset_train, gene_dataset_test, n_epochs=1)
