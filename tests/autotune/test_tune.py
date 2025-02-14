import os

import pytest
from ray import tune
from ray.tune import ResultGrid

from scvi import settings
from scvi.autotune import AutotuneExperiment, run_autotune
from scvi.data import synthetic_iid
from scvi.dataloaders import DataSplitter
from scvi.model import SCANVI, SCVI


def test_run_autotune_scvi_basic(save_path: str):
    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = run_autotune(
        SCVI,
        adata,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


def test_run_autotune_scvi_no_anndata(save_path: str, n_batches: int = 3):
    settings.logging_dir = save_path
    adata = synthetic_iid(n_batches=n_batches)
    SCVI.setup_anndata(adata, batch_key="batch")
    manager = SCVI._get_most_recent_anndata_manager(adata)

    datamodule = DataSplitter(manager)
    datamodule.n_vars = adata.n_vars
    datamodule.n_batch = n_batches

    experiment = run_autotune(
        SCVI,
        data=datamodule,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.parametrize("metric", ["Total", "Bio conservation", "iLISI"])
@pytest.mark.parametrize("model_cls", [SCVI, SCANVI])
def test_run_autotune_scvi_with_scib(model_cls, metric: str, save_path: str):
    # metric = "iLISI"
    # save_path = "."
    # Set an environment variable - only way it works in jax + ray
    os.environ["JAX_PLATFORMS"] = "cpu"
    os.environ["TUNE_DISABLE_STRICT_METRIC_CHECKING"] = "1"

    settings.logging_dir = save_path
    adata = synthetic_iid(batch_size=10, n_genes=10)
    if model_cls == SCANVI:
        model_cls.setup_anndata(
            adata,
            labels_key="labels",
            unlabeled_category="unknown",
            batch_key="batch",
        )
    else:
        model_cls.setup_anndata(
            adata,
            labels_key="labels",
            batch_key="batch",
        )

    experiment = run_autotune(
        model_cls,
        adata,
        metrics=[metric],
        mode="max",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


# def test_early_stopping():
#     # we use this temporarily to debug the scib-metrics callback
#     # (here we need to always allow extra metric in the vae)
#     n_epochs = 100
#
#     adata = synthetic_iid()
#     SCVI.setup_anndata(
#         adata,
#         batch_key="batch",
#         labels_key="labels",
#     )
#     model = SCVI(adata)
#     model.train(n_epochs, early_stopping=True, plan_kwargs={"lr": 0})
#     assert len(model.history["elbo_train"]) < n_epochs
