from random import randint

import scanpy as sc
from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler
from sklearn.metrics import silhouette_score

from scvi.autotune import Autotune
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_autotune():
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
    model, _ = tuner.run(metric="elbo_validation", scheduler=asha_scheduler)


def test_metric_function_dummy():
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
    model, analysis = tuner.run(metric="elbo_validation", scheduler=asha_scheduler)


def test_silhouette():
    def silhouette_metric(model):
        model.is_trained_ = True
        latent = model.get_latent_representation()
        model.is_trained_ = False
        adata.obsm["X_scVI"] = latent
        sc.pp.neighbors(adata, use_rep="X_scVI")
        sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.5)
        return silhouette_score(
            adata.obsm["X_scVI"],
            adata.obs["leiden_scVI"],
        )

    adata = synthetic_iid()

    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    metric_functions = {"silhouette_score": silhouette_metric}
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
    model, analysis = tuner.run(metric="elbo_validation", scheduler=asha_scheduler)
