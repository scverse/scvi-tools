from pytest import raises

from scvi import settings
from scvi.autotune import AutotuneExperiment
from scvi.autotune._experiment import _trainable
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_experiment_init(save_path: str):
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
    with raises(AttributeError):
        experiment.id = "new_id"

    assert hasattr(experiment, "data")
    assert experiment.data is not None
    assert experiment.data is adata
    with raises(AttributeError):
        experiment.data = "new_adata"

    assert hasattr(experiment, "setup_method_name")
    assert experiment.setup_method_name is not None
    assert experiment.setup_method_name == "setup_anndata"
    with raises(AttributeError):
        experiment.setup_method_name = "new_setup_method_name"

    assert hasattr(experiment, "setup_method_args")
    assert experiment.setup_method_args is not None
    with raises(AttributeError):
        experiment.setup_method_args = "new_setup_method_args"

    assert hasattr(experiment, "model_cls")
    assert experiment.model_cls is not None
    assert experiment.model_cls is SCVI
    with raises(AttributeError):
        experiment.model_cls = "new_model_cls"

    assert hasattr(experiment, "metrics")
    assert experiment.metrics is not None
    assert experiment.metrics == ["elbo_validation"]
    with raises(AttributeError):
        experiment.metrics = "new_metrics"

    assert hasattr(experiment, "mode")
    assert experiment.mode is not None
    assert experiment.mode == "min"
    with raises(AttributeError):
        experiment.mode = "new_mode"

    assert hasattr(experiment, "search_space")
    assert experiment.search_space is not None
    with raises(AttributeError):
        experiment.search_space = "new_search_space"

    assert hasattr(experiment, "num_samples")
    assert experiment.num_samples is not None
    assert experiment.num_samples == 1
    with raises(AttributeError):
        experiment.num_samples = 2

    assert hasattr(experiment, "scheduler")
    assert experiment.scheduler is not None
    with raises(AttributeError):
        experiment.scheduler = "new_scheduler"

    assert hasattr(experiment, "scheduler_kwargs")
    assert experiment.scheduler_kwargs is not None
    assert experiment.scheduler_kwargs == {}
    with raises(AttributeError):
        experiment.scheduler_kwargs = "new_scheduler_kwargs"

    assert hasattr(experiment, "searcher")
    assert experiment.searcher is not None
    with raises(AttributeError):
        experiment.searcher = "new_searcher"

    assert hasattr(experiment, "searcher_kwargs")
    assert experiment.searcher_kwargs is not None
    assert experiment.searcher_kwargs == {}
    with raises(AttributeError):
        experiment.searcher_kwargs = "new_searcher_kwargs"

    assert hasattr(experiment, "seed")
    assert experiment.seed is not None
    assert experiment.seed == settings.seed
    with raises(AttributeError):
        experiment.seed = 2

    assert hasattr(experiment, "resources")
    assert experiment.resources is not None
    assert experiment.resources == {}
    with raises(AttributeError):
        experiment.resources = "new_resources"

    assert hasattr(experiment, "name")
    assert experiment.name is not None
    assert experiment.name.startswith("scvi")
    with raises(AttributeError):
        experiment.name = "new_name"

    assert hasattr(experiment, "logging_dir")
    assert experiment.logging_dir is not None
    with raises(AttributeError):
        experiment.logging_dir = "new_logging_dir"

    with raises(AttributeError):
        _ = experiment.result_grid  # set after running the tuner


def test_experiment_no_setup_anndata():
    adata = synthetic_iid()

    with raises(ValueError):
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


def test_experiment_invalid_metrics():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(ValueError):
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
    with raises(ValueError):
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
    with raises(TypeError):
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


def test_experiment_invalid_mode():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(ValueError):
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
    with raises(ValueError):
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


def test_experiment_invalid_search_space():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(TypeError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space=None,
            num_samples=1,
        )
    with raises(ValueError):
        _ = AutotuneExperiment(
            SCVI,
            adata,
            metrics=["elbo_validation"],
            mode="min",
            search_space={},
            num_samples=1,
        )
    with raises(KeyError):
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


def test_experiment_invalid_num_samples():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(TypeError):
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


def test_experiment_invalid_scheduler():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(ValueError):
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
    with raises(TypeError):
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


def test_experiment_invalid_scheduler_kwargs():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(TypeError):
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


def test_experiment_invalid_searcher():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(ValueError):
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
    with raises(TypeError):
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


def test_experiment_invalid_searcher_kwargs():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(TypeError):
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


def test_experiment_invalid_seed():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    with raises(TypeError):
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


def test_experiment_get_tuner(save_path: str):
    from ray.tune import Tuner

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


def test_trainable(save_path: str):
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
