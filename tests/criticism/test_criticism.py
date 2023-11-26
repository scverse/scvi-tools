import pandas as pd
from anndata import AnnData
from sparse import GCXS
from xarray import Dataset

from scvi.criticism import PosteriorPredictiveCheck as PPC
from scvi.data import synthetic_iid
from scvi.model import SCVI


def get_ppc_with_samples(adata: AnnData, n_samples: int = 2):
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
    ppc = PPC(adata, models_dict, n_samples=n_samples)
    return ppc, models_dict


def test_ppc_init():
    adata = synthetic_iid()
    ppc, models_dict = get_ppc_with_samples(adata, n_samples=42)
    assert isinstance(ppc.raw_counts, GCXS)
    assert isinstance(ppc.samples_dataset, Dataset)
    assert ppc.n_samples == 42
    assert ppc.models is models_dict
    assert ppc.metrics == {}
    assert ppc.samples_dataset.model1.shape == (400, 100, 42)
    assert ppc.samples_dataset.model2.shape == (400, 100, 42)


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


def test_ppc_de(n_genes: int = 200):
    adata = synthetic_iid(n_genes=n_genes)
    ppc, _ = get_ppc_with_samples(adata, n_samples=4)

    # Use a high thresh for simulated data
    ppc.differential_expression(de_groupby="labels", p_val_thresh=0.7)
