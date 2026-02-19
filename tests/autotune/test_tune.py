import pytest

from scvi import settings
from scvi.data import synthetic_iid
from scvi.dataloaders import DataSplitter
from scvi.external import TOTALANVI
from scvi.model import MULTIVI, SCANVI, SCVI, TOTALVI


@pytest.mark.autotune
@pytest.mark.parametrize("save_checkpoints", [True, False])
@pytest.mark.parametrize("metric", ["elbo_validation"])
def test_run_autotune_scvi_basic_adata(save_checkpoints: bool, metric: str, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = run_autotune(
        SCVI,
        adata,
        metrics=[metric],
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
        save_checkpoints=save_checkpoints,
        scheduler="asha",
        searcher="hyperopt",
        ignore_reinit_error=True,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.autotune
@pytest.mark.parametrize("save_checkpoints", [True, False])
def test_run_autotune_scvi_basic_mdata(save_checkpoints: bool, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

    settings.logging_dir = save_path
    mdata = synthetic_iid(return_mudata=True)
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={
            "rna_layer": "rna",
            "atac_layer": "accessibility",
            "protein_layer": "protein_expression",
        },
    )

    experiment = run_autotune(
        MULTIVI,
        mdata,
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
        save_checkpoints=save_checkpoints,
        scheduler="asha",
        searcher="hyperopt",
        ignore_reinit_error=True,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.autotune
@pytest.mark.parametrize("n_batches", [3])
def test_run_autotune_scvi_no_anndata(n_batches: int, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

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
        ignore_reinit_error=True,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.autotune
@pytest.mark.parametrize("metric", ["Total"])
@pytest.mark.parametrize("model_cls", [SCVI, TOTALVI])
@pytest.mark.parametrize("solver", ["randomized"])
def test_run_autotune_scvi_with_scib_adata(model_cls, metric: str, solver: str, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

    settings.logging_dir = save_path
    adata = synthetic_iid()
    if model_cls == SCANVI:
        model_cls.setup_anndata(
            adata,
            labels_key="labels",
            unlabeled_category="unknown",
            batch_key="batch",
        )
    elif model_cls == SCVI:
        model_cls.setup_anndata(
            adata,
            labels_key="labels",
            batch_key="batch",
        )
    elif model_cls == TOTALVI:
        model_cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
        )
    elif model_cls == TOTALANVI:
        model_cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
            labels_key="labels",
            unlabeled_category="label_0",
        )
    else:
        raise ValueError("No Model")

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
        scib_subsample_rows=100,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
        local_mode=True,
        ignore_reinit_error=True,
        solver=solver,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.autotune
@pytest.mark.parametrize("metric", ["Total"])
@pytest.mark.parametrize("model_cls", [MULTIVI, TOTALVI])
@pytest.mark.parametrize("solver", ["randomized"])
def test_run_autotune_scvi_with_scib_mdata(model_cls, metric: str, solver: str, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

    settings.logging_dir = save_path
    mdata = synthetic_iid(return_mudata=True)
    if model_cls == MULTIVI:
        model_cls.setup_mudata(
            mdata,
            batch_key="batch",
            modalities={
                "rna_layer": "rna",
                "atac_layer": "accessibility",
                "protein_layer": "protein_expression",
            },
        )
    elif model_cls == TOTALVI:
        model_cls.setup_mudata(
            mdata,
            batch_key="batch",
            modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
        )
    else:
        raise ValueError("No Model")

    experiment = run_autotune(
        model_cls,
        mdata,
        metrics=[metric],
        mode="max",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([5, 10]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        scib_subsample_rows=100,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
        local_mode=True,
        ignore_reinit_error=True,
        solver=solver,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.autotune
@pytest.mark.parametrize("metric", ["iLISI"])
def test_run_autotune_scvi_with_scib_ext_indices(metric: str, save_path: str):
    from ray import tune
    from ray.tune import ResultGrid

    from scvi.autotune import AutotuneExperiment, run_autotune

    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCANVI.setup_anndata(
        adata,
        labels_key="labels",
        unlabeled_category="unknown",
        batch_key="batch",
    )

    experiment = run_autotune(
        SCANVI,
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
        scib_indices_list=[1, 2, 3],
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
        local_mode=True,
        ignore_reinit_error=True,
        solver="randomized",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)
