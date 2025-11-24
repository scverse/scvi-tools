from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from sparse import GCXS
from xarray import Dataset

from scvi.criticism import PosteriorPredictiveCheck as PPC
from scvi.criticism import create_criticism_report
from scvi.data import synthetic_iid
from scvi.model import SCVI

if TYPE_CHECKING:
    from anndata import AnnData


def get_ppc_with_samples(adata: AnnData, n_samples: int = 2, indices: list[int] | None = None):
    # create and train models
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model1 = SCVI(adata, n_latent=5)
    model1.train(1)

    bdata = adata.copy()
    SCVI.setup_anndata(
        bdata,
        batch_key="batch",
    )
    model2 = SCVI(bdata, n_latent=5)
    model2.train(1)

    models_dict = {"model1": model1, "model2": model2}
    ppc = PPC(adata, models_dict, n_samples=n_samples, indices=indices)
    return ppc, models_dict


@pytest.mark.optional
def test_ppc_init(save_path):
    adata = synthetic_iid()
    ppc, models_dict = get_ppc_with_samples(adata, n_samples=42)
    model1 = models_dict["model1"]
    model_path1 = os.path.join(save_path, "model1")
    model1.save(model_path1, save_anndata=True, overwrite=True)
    model2 = models_dict["model2"]
    model_path2 = os.path.join(save_path, "model2")
    model2.save(model_path2, save_anndata=True, overwrite=True)
    create_criticism_report(model1, save_folder=model_path1)
    create_criticism_report(model2, save_folder=model_path2)
    assert isinstance(ppc.raw_counts, GCXS)
    assert isinstance(ppc.samples_dataset, Dataset)
    assert ppc.n_samples == 42
    assert ppc.models is models_dict
    assert ppc.metrics == {}
    assert ppc.samples_dataset.model1.shape == (adata.n_obs, 100, 42)
    assert ppc.samples_dataset.model2.shape == (adata.n_obs, 100, 42)

    ppc, models_dict = get_ppc_with_samples(adata, n_samples=42, indices=np.arange(100))
    assert isinstance(ppc.raw_counts, GCXS)
    assert isinstance(ppc.samples_dataset, Dataset)
    assert ppc.n_samples == 42
    assert ppc.models is models_dict
    assert ppc.metrics == {}
    assert ppc.samples_dataset.model1.shape == (100, 100, 42)
    assert ppc.samples_dataset.model2.shape == (100, 100, 42)


def test_ppc_cv():
    adata = synthetic_iid(n_genes=10)
    ppc, _ = get_ppc_with_samples(adata)

    ppc.coefficient_of_variation("cells")
    ppc.coefficient_of_variation("features")

    assert list(ppc.metrics.keys()) == ["cv_gene", "cv_cell"]

    assert isinstance(ppc.metrics["cv_cell"], pd.DataFrame)
    assert ppc.metrics["cv_cell"].columns.tolist() == ["model1", "model2", "Raw"]
    assert ppc.metrics["cv_cell"].index.equals(adata.obs_names)

    assert isinstance(ppc.metrics["cv_gene"], pd.DataFrame)
    assert ppc.metrics["cv_gene"].columns.tolist() == ["model1", "model2", "Raw"]
    assert ppc.metrics["cv_gene"].index.equals(adata.var_names)


def test_ppc_calibration():
    adata = synthetic_iid(n_genes=10)
    ppc, _ = get_ppc_with_samples(adata, n_samples=4)

    ppc.calibration_error()

    assert list(ppc.metrics.keys()) == ["calibration"]

    assert isinstance(ppc.metrics["calibration"], pd.DataFrame)
    assert ppc.metrics["calibration"].columns.tolist() == ["model1", "model2"]


def test_ppc_zero_fraction():
    adata = synthetic_iid(n_genes=10)
    ppc, _ = get_ppc_with_samples(adata)

    ppc.zero_fraction()

    assert list(ppc.metrics.keys()) == ["zero_fraction"]

    assert isinstance(ppc.metrics["zero_fraction"], pd.DataFrame)
    assert ppc.metrics["zero_fraction"].columns.tolist() == ["model1", "model2", "Raw"]


@pytest.mark.parametrize("n_genes", [200])
def test_ppc_de(n_genes: int):
    adata = synthetic_iid(n_genes=n_genes)
    ppc, _ = get_ppc_with_samples(adata, n_samples=4)

    # Use a high thresh for simulated data
    ppc.differential_expression(de_groupby="labels", p_val_thresh=0.7)
