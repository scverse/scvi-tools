from functools import partial

import numpy as np
import pandas as pd
from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune, silhouette_metric_labels_batch, tune_scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_tune_scvi(save_path):
    adata = synthetic_iid()
    tune_scvi(adata, 2, run_kwargs=dict(local_dir=save_path))


def test_tune_scvi_train(save_path):
    adata = synthetic_iid()
    tune_scvi(adata, 2, run_kwargs=dict(local_dir=save_path, train=True))


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
    print(model)


def test_metric_function_dummy(save_path):
    adata = synthetic_iid()
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    metric_functions = {"dummy": lambda x: np.random.randint(1, 10)}
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


def test_both_covariates(save_path):
    adata = synthetic_iid()
    adata.obs["continious_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

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
    continious_covariates = ["continious_cov"]
    categorical_covariates = ["categorical_cov"]
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
        continious_covariates=continious_covariates,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )


def test_effect_covariate(save_path):
    adata = synthetic_iid()
    adata.obs["continious_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

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
    categorical_covariates = ["categorical_cov"]
    continious_covariates = ["continious_cov"]
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
        continious_covariates=continious_covariates,
        test_effect_covariates=False,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )


def test_single_covariate(save_path):
    adata = synthetic_iid()
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

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
    categorical_covariates = ["categorical_cov"]
    num_epochs = 2
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )


def test_hvg_filter(save_path):
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
        test_effect_hvg=True,
        top_hvg=[1000, 2000],
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )


def test_hvg_no_filter(save_path):
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
        top_hvg=[1000, 2000],
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )


def test_complex(save_path):
    adata = synthetic_iid()
    adata.obs["continious_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )
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
    continious_covariates = ["continious_cov"]
    categorical_covariates = ["categorical_cov", "batch"]
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        test_effect_hvg=True,
        top_hvg=[2000, 5000],
        continious_covariates=continious_covariates,
        categorical_covariates=categorical_covariates,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
