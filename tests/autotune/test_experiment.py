import pytest

from scvi import settings
from scvi.data import synthetic_iid
from scvi.model import SCVI


@pytest.mark.autotune
def test_experiment_init(save_path: str):
    from scvi.autotune import AutotuneExperiment

    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = AutotuneExperiment(
        SCVI,
        adata,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": [1, 2],
            }
        },
        num_samples=1,
    )
    assert hasattr(experiment, "id")
    assert experiment.id is not None
    assert isinstance(experiment.id, str)
    with pytest.raises(AttributeError):
        experiment.id = "new_id"

    assert hasattr(experiment, "data")
    assert experiment.data is not None
    assert experiment.data is adata
    with pytest.raises(AttributeError):
        experiment.data = "new_adata"

    assert hasattr(experiment, "setup_method_name")
    assert experiment.setup_method_name is not None
    assert experiment.setup_method_name == "setup_anndata"
    with pytest.raises(AttributeError):
        experiment.setup_method_name = "new_setup_method_name"

    assert hasattr(experiment, "setup_method_args")
    assert experiment.setup_method_args is not None
    with pytest.raises(AttributeError):
        experiment.setup_method_args = "new_setup_method_args"

    assert hasattr(experiment, "model_cls")
    assert experiment.model_cls is not None
    assert experiment.model_cls is SCVI
    with pytest.raises(AttributeError):
        experiment.model_cls = "new_model_cls"

    assert hasattr(experiment, "metrics")
    assert experiment.metrics is not None
    assert experiment.metrics == ["elbo_validation"]
    with pytest.raises(AttributeError):
        experiment.metrics = "new_metrics"

    assert hasattr(experiment, "mode")
    assert experiment.mode is not None
    assert experiment.mode == "min"
    with pytest.raises(AttributeError):
        experiment.mode = "new_mode"

    assert hasattr(experiment, "search_space")
    assert experiment.search_space is not None
    with pytest.raises(AttributeError):
        experiment.search_space = "new_search_space"

    assert hasattr(experiment, "num_samples")
    assert experiment.num_samples is not None
    assert experiment.num_samples == 1
    with pytest.raises(AttributeError):
        experiment.num_samples = 2

    assert hasattr(experiment, "scheduler")
    assert experiment.scheduler is not None
    with pytest.raises(AttributeError):
        experiment.scheduler = "new_scheduler"

    assert hasattr(experiment, "scheduler_kwargs")
    assert experiment.scheduler_kwargs is not None
    assert experiment.scheduler_kwargs == {}
    with pytest.raises(AttributeError):
        experiment.scheduler_kwargs = "new_scheduler_kwargs"

    assert hasattr(experiment, "searcher")
    assert experiment.searcher is not None
    with pytest.raises(AttributeError):
        experiment.searcher = "new_searcher"

    assert hasattr(experiment, "searcher_kwargs")
    assert experiment.searcher_kwargs is not None
    assert experiment.searcher_kwargs == {}
    with pytest.raises(AttributeError):
        experiment.searcher_kwargs = "new_searcher_kwargs"

    assert hasattr(experiment, "seed")
    assert experiment.seed is not None
    assert experiment.seed == settings.seed
    with pytest.raises(AttributeError):
        experiment.seed = 2

    assert hasattr(experiment, "resources")
    assert experiment.resources is not None
    assert experiment.resources == {}
    with pytest.raises(AttributeError):
        experiment.resources = "new_resources"

    assert hasattr(experiment, "name")
    assert experiment.name is not None
    assert experiment.name.startswith("scvi")
    with pytest.raises(AttributeError):
        experiment.name = "new_name"

    assert hasattr(experiment, "logging_dir")
    assert experiment.logging_dir is not None
    with pytest.raises(AttributeError):
        experiment.logging_dir = "new_logging_dir"

    with pytest.raises(AttributeError):
        _ = experiment.result_grid  # set after running the tuner


@pytest.mark.autotune
def test_experiment_no_setup_anndata():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()

    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )


@pytest.mark.autotune
def test_experiment_invalid_metrics():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=[],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )
    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=None,
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )
    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics={"elbo_validation": None},
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )


@pytest.mark.autotune
def test_experiment_invalid_mode():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="minimum",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )
    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode=None,
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )


@pytest.mark.autotune
def test_experiment_invalid_search_space():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space=None,
            num_samples=1,
        )
    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={},
            num_samples=1,
        )
    with pytest.raises(KeyError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "invalid_key_here": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
        )


@pytest.mark.autotune
def test_experiment_invalid_num_samples():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=None,
        )


@pytest.mark.autotune
def test_experiment_invalid_scheduler():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            scheduler="invalid option",
        )
    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            scheduler=[],
        )


@pytest.mark.autotune
def test_experiment_invalid_scheduler_kwargs():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            scheduler_kwargs="invalid type",
        )


@pytest.mark.autotune
def test_experiment_invalid_searcher():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            searcher="invalid option",
        )
    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            searcher=[],
        )


@pytest.mark.autotune
def test_experiment_invalid_searcher_kwargs():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            searcher_kwargs="invalid type",
        )


@pytest.mark.autotune
def test_experiment_invalid_seed():
    from scvi.autotune import AutotuneExperiment

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with pytest.raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={
                "model_params": {
                    "n_hidden": [1, 2],
                }
            },
            num_samples=1,
            seed="invalid type",
        )


@pytest.mark.autotune
def test_experiment_get_tuner(save_path: str):
    from ray.tune import Tuner

    from scvi.autotune import AutotuneExperiment

    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = AutotuneExperiment(
        SCVI,
        adata,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": [1, 2],
            }
        },
        num_samples=1,
    )
    tuner = experiment.get_tuner()
    assert isinstance(tuner, Tuner)


@pytest.mark.autotune
def test_trainable(save_path: str):
    from scvi.autotune import AutotuneExperiment
    from scvi.autotune._experiment import _trainable

    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = AutotuneExperiment(
        SCVI,
        adata,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": [1, 2],
            }
        },
        num_samples=1,
    )
    sample = {
        "model_params": {
            "n_hidden": 1,
        }
    }
    _trainable(sample, experiment)
