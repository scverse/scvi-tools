from functools import partial
from random import randint

from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune, silhouette_metric_labels_batch, tune_scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_tune_scvi(save_path):
    adata = synthetic_iid()
    tune_scvi(adata, 2, run_kwargs=dict(local_dir=save_path))


def test_autotune(save_path):
    adata = synthetic_iid()
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    model_config = {"dropout_rate": loguniform(1e-4, 1e-1)}
    plan_config = {"lr": loguniform(1e-4, 1e-1)}
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="elbo_validation", scheduler=asha_scheduler, local_dir=save_path
    )


def test_metric_function_dummy(save_path):
    adata = synthetic_iid()
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    metric_functions = {"dummy": lambda x: randint(1, 10)}
    model_config = {"dropout_rate": loguniform(1e-4, 1e-1)}
    plan_config = {"lr": loguniform(1e-4, 1e-1)}
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="elbo_validation", scheduler=asha_scheduler, local_dir=save_path
    )


def test_silhouette(save_path):

    adata = synthetic_iid()

    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    metric_functions = {
        "silhouette_score": partial(
            silhouette_metric_labels_batch, labels_key="labels", batch_key="batch"
        )
    }
    model_config = {"dropout_rate": loguniform(1e-4, 1e-1)}
    plan_config = {"lr": loguniform(1e-4, 1e-1)}
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
