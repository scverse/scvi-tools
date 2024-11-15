from ray import tune
from ray.tune import ResultGrid

from scvi import settings
from scvi.autotune import AutotuneExperiment, run_autotune
from scvi.data import synthetic_iid
from scvi.model import SCVI


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
        num_samples=1,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


def test_run_autotune_scvi_no_anndata(save_path: str, n_batches: int = 3):
    from scvi.dataloaders import DataSplitter

    settings.logging_dir = save_path
    adata = synthetic_iid(n_batches=n_batches)
    SCVI.setup_anndata(adata, batch_key="batch")
    manager = SCVI._get_most_recent_anndata_manager(adata)

    datamodule = DataSplitter(manager)
    datamodule.n_vars = adata.n_vars
    datamodule.n_batch = n_batches

    experiment = run_autotune(
        SCVI,
        datamodule,
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
        num_samples=1,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)
