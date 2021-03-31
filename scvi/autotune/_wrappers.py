from ray.tune import loguniform
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune
from scvi.model import SCVI


def tune_scvi(adata):
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
    return tuner.run(metric="elbo_validation", scheduler=asha_scheduler)
