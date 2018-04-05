#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

from run_benchmarks import run_benchmarks
from scvi.dataset import load_dataset


def test_benchmark():
    run_benchmarks(load_dataset("synthetic"), n_epochs=1)
