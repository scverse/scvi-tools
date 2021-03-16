from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler, PopulationBasedTraining

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
    tuner.run(metric="elbo_validation", scheduler=asha_scheduler)


def test_pbt():
    adata = synthetic_iid()
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    model_config = {"dropout_rate": 1e-3}
    plan_config = {"lr": 1e-3}
    tuner = Autotune(
        adata,
        SCVI,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
    )
    pbt_scheduler = PopulationBasedTraining(
        perturbation_interval=4,
        hyperparam_mutations={
            "lr": loguniform(1e-4, 1e-1),
            "dropout_rate": loguniform(1e-4, 1e-1),
        },
    )
    tuner.run(metric="elbo_validation", scheduler=pbt_scheduler)


def test_metric_functions():
    pass
