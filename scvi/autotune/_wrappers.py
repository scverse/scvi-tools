import anndata
from ray.tune import choice, loguniform
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune
from scvi.model import SCVI


def tune_scvi(adata: anndata.AnnData, n_epochs, run_kwargs):
    """
    Tune scvi with defaults for `tune.run` and return the best model.

    Parameters
    ----------
    adata
        AnnData object we will tune the model on.
    n_epochs
        Max epochs for `tune.run`.

    Returns
    -------
    A tuple with the best model object and tune `Analysis` object
    """
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    metric_functions = {}
    model_config = {
        "dropout_rate": loguniform(1e-4, 1e-1),
        "n_layers": choice([i for i in range(1, 5)]),
        "n_hidden": choice([64, 128, 256]),
        "n_latent": choice([i for i in range(5, 15)]),
    }
    plan_config = {"lr": loguniform(1e-4, 1e-1)}
    tuner = Autotune(
        adata,
        SCVI,
        metric_functions=metric_functions,
        training_metrics=metrics,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        num_epochs=n_epochs,
    )
    asha_scheduler = ASHAScheduler(max_t=n_epochs, grace_period=1, reduction_factor=2)
    return tuner.run(metric="elbo_validation", scheduler=asha_scheduler, **run_kwargs)
