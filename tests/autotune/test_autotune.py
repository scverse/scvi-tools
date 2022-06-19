from functools import partial

import numpy as np
import pandas as pd
from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune, silhouette_metric_labels_batch, tune_scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI

MODEL_CONFIG = {"dropout_rate": loguniform(1e-4, 1e-1)}
PLAN_CONFIG = {"lr": loguniform(1e-4, 1e-1)}
NUM_EPOCHS = 2
METRICS = ["elbo_validation", "reconstruction_loss_validation"]
METRIC_FUNCTIONS = {
    "silhouette_score": partial(
        silhouette_metric_labels_batch, labels_key="labels", batch_key="batch"
    )
}


def test_tune_scvi(save_path):
    adata = synthetic_iid()
    tune_scvi(adata, NUM_EPOCHS, run_kwargs=dict(local_dir=save_path))


def test_tune_scvi_train(save_path):
    adata = synthetic_iid()
    tune_scvi(adata, NUM_EPOCHS, run_kwargs=dict(local_dir=save_path, train=True))


def test_autotune(save_path):
    adata = synthetic_iid()
    metrics = METRICS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
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


def test_autotune_setupargs(save_path):
    adata = synthetic_iid()
    metrics = METRICS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    setup_args = {"batch_key": "batch", "categorical_covariate_keys": ["labels"]}
    num_epochs = NUM_EPOCHS
    tuner = Autotune(
        adata,
        SCVI,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        setup_args=setup_args,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="elbo_validation", scheduler=asha_scheduler, local_dir=save_path
    )
    print(model)


def test_metric_function_dummy(save_path):
    adata = synthetic_iid()
    metric_functions = {"dummy": lambda x: np.random.randint(1, 10)}
    metrics = METRICS
    metric_functions = {"dummy": lambda x: np.random.randint(1, 10)}
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
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
    print(model)


def test_silhouette(save_path):

    adata = synthetic_iid()

    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
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
    print(model)


def test_both_covariates(save_path):
    adata = synthetic_iid()
    adata.obs["continuous_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    continuous_covariates = ["continuous_cov"]
    categorical_covariates = ["categorical_cov"]
    num_epochs = NUM_EPOCHS
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
        continuous_covariates=continuous_covariates,
        test_effect_covariates=True,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
    print(model)


def test_covariates_elbo(save_path):
    adata = synthetic_iid()
    adata.obs["continuous_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    continuous_covariates = ["continuous_cov"]
    categorical_covariates = ["categorical_cov"]
    num_epochs = NUM_EPOCHS
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
        continuous_covariates=continuous_covariates,
        test_effect_covariates=True,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="elbo_validation",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
    print(model)


def test_effect_covariate(save_path):
    adata = synthetic_iid()
    adata.obs["continuous_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    categorical_covariates = ["categorical_cov"]
    continuous_covariates = ["continuous_cov"]
    num_epochs = NUM_EPOCHS
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        categorical_covariates=categorical_covariates,
        continuous_covariates=continuous_covariates,
        test_effect_covariates=True,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
    print(model)


def test_single_covariate(save_path):
    adata = synthetic_iid()
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )

    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    categorical_covariates = ["categorical_cov"]
    num_epochs = NUM_EPOCHS
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
    print(model)


def test_hvg_filter(save_path):
    adata = synthetic_iid()
    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
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
    print(model)


def test_hvg_no_filter(save_path):
    adata = synthetic_iid()
    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
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
    print(model)


def test_complex(save_path):
    adata = synthetic_iid()
    adata.obs["continuous_cov"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["categorical_cov"] = pd.Categorical(
        ["cov_1", "cov_2"] * (round(adata.n_obs / 2))
    )
    metrics = METRICS
    metric_functions = METRIC_FUNCTIONS
    model_config = MODEL_CONFIG
    plan_config = PLAN_CONFIG
    num_epochs = NUM_EPOCHS
    continuous_covariates = ["continuous_cov"]
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
        continuous_covariates=continuous_covariates,
        categorical_covariates=categorical_covariates,
        test_effect_covariates=True,
    )
    asha_scheduler = ASHAScheduler(max_t=num_epochs, grace_period=1, reduction_factor=2)
    model, analysis = tuner.run(
        metric="silhouette_score",
        mode="max",
        scheduler=asha_scheduler,
        local_dir=save_path,
    )
    print(model)
