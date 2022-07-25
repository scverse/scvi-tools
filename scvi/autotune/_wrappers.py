from typing import Optional, Union

import anndata
from ray import tune
from ray.tune import ExperimentAnalysis
from ray.tune.schedulers import ASHAScheduler

from scvi.autotune import Autotune
from scvi.model import SCVI
from scvi.model.base import BaseModelClass


def tune_scvi(
    adata: anndata.AnnData,
    n_epochs: int,
    run_kwargs: Optional[dict] = None,
) -> Union[BaseModelClass, ExperimentAnalysis]:
    """
    Tune scvi with defaults for `tune.run` and return the best model.

    The following hyperparameters will be compared:
    - dropout_rate: `loguniform(1e-4, 1e-1)`
    - n_layers: `random integer between 1 and 5`
    - n_hidden: `random integer 64, 128 and 256`
    - n_latent: `random integer between 20 and 50 in increments of 5`

    The default metrics used to optimize the model are `elbo_validation` and `reconstruction_loss_validation`.

    Parameters
    ----------
    adata
        AnnData object we will tune the model on.
    n_epochs
        Max epochs for `tune.run`.
    run_kwargs
        Extra keyword arguments for `tune.run`.

    Returns
    -------
    A tuple with the best model object and tune `Analysis` object

    Notes
    -------
    This is a wrapper with limited functionality. If you want to have more control over the hyperparameters,
    their values and the metrics used to optimize the model, please refer to :class:`~scvi.autotune.Autotune`
    """
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
    ]
    model_config = {
        "dropout_rate": tune.loguniform(1e-4, 1e-1),
        "n_layers": tune.choice([i for i in range(1, 5)]),
        "n_hidden": tune.choice([64, 128, 256]),
        "n_latent": tune.choice(list(range(20, 55, 5))),
    }
    plan_config = {"lr": tune.loguniform(1e-4, 1e-1)}
    tuner = Autotune(
        adata,
        SCVI,
        model_hyperparams=model_config,
        plan_hyperparams=plan_config,
        num_epochs=n_epochs,
        training_metrics=metrics,
    )
    asha_scheduler = ASHAScheduler(max_t=n_epochs, grace_period=1, reduction_factor=2)
    return tuner.run(metric="elbo_validation", scheduler=asha_scheduler, **run_kwargs)
