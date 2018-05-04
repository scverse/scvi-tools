#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""
from run_benchmarks import run_benchmarks
from scvi.benchmark import run_benchmarks_classification
from scvi.models import VAEC, VAE, SVAEC


def test_synthetic_1():
    run_benchmarks("synthetic", n_epochs=1, use_batches=True, model=VAE)
    run_benchmarks_classification("synthetic", n_epochs=1)


def test_synthetic_2():
    run_benchmarks("synthetic", n_epochs=1, model=SVAEC, benchmark=True)


def test_cortex():
    run_benchmarks("cortex", n_epochs=1, model=VAEC)


def test_brain_large():
    run_benchmarks("brain_large", n_epochs=1, use_batches=False, tt_split=0.5, unit_test=True)


def test_retina():
    run_benchmarks("retina", n_epochs=1, show_batch_mixing=False, unit_test=True)
