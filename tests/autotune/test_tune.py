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
            "model_args": {
                "n_hidden": tune.choice([1, 2]),
            }
        },
        max_epochs=1,
        num_samples=1,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)
