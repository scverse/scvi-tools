#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""
from scvi.benchmark import run_benchmarks
from scvi.benchmark import run_benchmarks_classification
from scvi.models import VAEC, VAE, SVAEC


def test_synthetic_1():
    run_benchmarks("synthetic", n_epochs=1, use_batches=True, model=VAE)
    run_benchmarks_classification("synthetic", n_epochs=1, n_epochs_classifier=1)


def test_synthetic_2():
    run_benchmarks("synthetic", n_epochs=1, model=SVAEC, benchmark=True)


def test_cortex():
    run_benchmarks("cortex", n_epochs=1, model=VAEC)


def test_brain_large():
    run_benchmarks("brain_large", n_epochs=1, use_batches=False, tt_split=0.5, save_path='tests/data/')


def test_retina():
    run_benchmarks("retina", n_epochs=1, show_batch_mixing=False, save_path='tests/data/')


def test_cbmc():
    run_benchmarks("cbmc", n_epochs=1, show_batch_mixing=False, save_path='tests/data/')


def test_brain_small():
    run_benchmarks("brain_small", n_epochs=1, show_batch_mixing=False, save_path='tests/data/')


def test_hemato():
    run_benchmarks("hemato", n_epochs=1, show_batch_mixing=False, save_path='tests/data/HEMATO/')


def test_pbmc():
    run_benchmarks("pbmc", n_epochs=1, show_batch_mixing=False, save_path='tests/data/PBMC/')


def test_loom():
    run_benchmarks(dataset_name="retina.loom", n_epochs=1, show_batch_mixing=False, save_path='tests/data/')
