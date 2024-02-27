from __future__ import annotations

import logging
from typing import Any, Literal

from ray.tune import Tuner

from scvi._types import AnnOrMuData
from scvi.autotune._experiment import AutotuneExperiment
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


def _trainable(
    search_space: dict[str, callable],
    experiment: AutotuneExperiment,
) -> None:
    """Implements a Ray Tune trainable function for an :class:`~scvi.autotune.AutotuneExperiment`.

    Parameters
    ----------
    search_space
        Hyperparameter configuration to evaluate. Note: this is different from
        :attr:`~scvi.autotune.AutotuneExperiment.search_space` in that it is a single configuration
        sampled from the search space, not the specification of the search space itself.
    experiment
        :class:`~scvi.autotune.AutotuneExperiment` to evaluate.

    Notes
    -----
    See the Ray Tune
    `documentation <https://docs.ray.io/en/latest/tune/api/trainable.html#function-trainable-api>`_
    for more details.
    """
    from os.path import join

    from lightning.pytorch.callbacks import Callback
    from lightning.pytorch.loggers import TensorBoardLogger
    from ray.train import get_context
    from ray.tune.integration.pytorch_lightning import TuneReportCheckpointCallback

    from scvi import settings

    # This is to get around lightning import changes
    tune_callback_cls = type(
        "_TuneReportCheckpointCallback",
        (TuneReportCheckpointCallback, Callback),
        {},
    )

    callbacks = [
        tune_callback_cls(metrics=experiment.metrics, on="validation_end", save_checkpoints=False)
    ]
    log_dir = join(
        experiment.logging_dir,
        experiment.name,
        f"{get_context().get_trial_name()}_tensorboard",
    )

    model_args = search_space.get("model_args", {})
    train_args = search_space.get("train_args", {})
    train_args = {
        "max_epochs": experiment.max_epochs,
        "accelerator": "auto",
        "devices": "auto",
        "check_val_every_n_epoch": 1,
        "enable_progress_bar": False,
        "logger": TensorBoardLogger(log_dir),
        "callbacks": callbacks,
        **train_args,
    }

    settings.seed = experiment.seed
    getattr(experiment.model_cls, experiment.setup_method_name)(
        experiment.adata,
        **experiment.setup_method_args,
    )
    model = experiment.model_cls(experiment.adata, **model_args)
    model.train(**train_args)


def _configure_tuner(experiment: AutotuneExperiment) -> Tuner:
    """Configure a :class:`~ray.tune.Tuner` for an :class:`~scvi.autotune.AutotuneExperiment`."""
    from ray.train import RunConfig
    from ray.tune import with_parameters, with_resources
    from ray.tune.tune_config import TuneConfig

    trainable = with_parameters(_trainable, experiment=experiment)
    trainable = with_resources(trainable, resources=experiment.resources)

    tune_config = TuneConfig(
        scheduler=experiment.scheduler,
        search_alg=experiment.searcher,
        num_samples=experiment.num_samples,
    )
    run_config = RunConfig(
        name=experiment.name,
        storage_path=experiment.logging_dir,
        local_dir=experiment.logging_dir,
        log_to_file=True,
        verbose=1,
    )
    return Tuner(
        trainable=trainable,
        param_space=experiment.search_space,
        tune_config=tune_config,
        run_config=run_config,
    )


def run_autotune(
    model_cls: BaseModelClass,
    adata: AnnOrMuData,
    metrics: str | list[str],
    mode: Literal["min", "max"],
    search_space: dict[str, dict[Literal["model_args", "train_args"], dict[str, Any]]],
    num_samples: int,
    max_epochs: int | None = None,
    scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
    searcher: Literal["hyperopt", "random"] = "hyperopt",
    seed: int | None = None,
    resources: dict[Literal["cpu", "gpu", "memory"], float] | None = None,
    experiment_name: str | None = None,
    logging_dir: str | None = None,
    scheduler_kwargs: dict | None = None,
    searcher_kwargs: dict | None = None,
) -> AutotuneExperiment:
    """``BETA`` Run a hyperparameter sweep.

    Parameters
    ----------
    model_cls
        Model class on which to tune hyperparameters.
    adata
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup with ``model_cls``.
    metrics
        Either a single metric or a list of metrics to track during the experiment. If a list is
        provided, the primary metric will be the first element in the list.
    mode
        Optimization mode for the primary metric. One of ``"min"`` or ``"max"``.
    search_space
        Dictionary of hyperparameter names and their respective search spaces. See
        the `API <https://docs.ray.io/en/latest/tune/api/search_space.html>`_ for available
        search specifications. Must only contain the following top-level keys:

        * ``"model_args"``: parameters to pass to the model constructor.
        * ``"train_args"``: parameters to pass to the model's ``train`` method.
    num_samples
        Total number of hyperparameter configurations to sample from the search space. Passed
        into :class:`~ray.tune.tune_config.TuneConfig`.
    max_epochs
        Maximum number of epochs to train hyperparameter configurations.
    scheduler
        Ray Tune scheduler to use. One of the following:

        * ``"asha"``: :class:`~ray.tune.schedulers.AsyncHyperBandScheduler`
        * ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
        * ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
        * ``"pbt"``: :class:`~ray.tune.schedulers.PopulationBasedTraining`
        * ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`

        Note that that not all schedulers are compatible with all search algorithms. See the
        `documentation <https://docs.ray.io/en/latest/tune/key-concepts.html#schedulers>`_
        for more details.
    searcher
        Ray Tune search algorithm to use. One of the following:

        * ``"hyperopt"``: :class:`~ray.tune.search.hyperopt.HyperOptSearch`
        * ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`
    seed
        Random seed to use for the experiment. Propagated to :attr:`~scvi.settings.seed` and
        search algorithms. If not provided, defaults to :attr:`~scvi.settings.seed`.
    resources
        Dictionary of resources to allocate per trial in the experiment. Available keys
        include:

        * ``"cpu"``: number of CPU cores
        * ``"gpu"``: number of GPUs
        * ``"memory"``: amount of memory

        If not provided, defaults to using all available resources. Note that fractional
        allocations are supported.
    experiment_name
        Name of the experiment, used for logging purposes. Defaults to a unique ID.
    logging_dir
        Base directory to store experiment logs. Defaults to :attr:``scvi.settings.logging_dir``.
    scheduler_kwargs
        Additional keyword arguments to pass to the scheduler.
    searcher_kwargs
        Additional keyword arguments to pass to the search algorithm.

    Returns
    -------
    :class:`~scvi.autotune.AutotuneExperiment` object containing the results of the hyperparameter
    sweep in :attr:`~scvi.autotune.AutotuneExperiment.result_grid`.

    Notes
    -----
    Lifecycle: beta

    See Also
    --------
    :class:`~scvi.autotune.AutotuneExperiment`
    """
    from ray import init

    experiment = AutotuneExperiment(
        model_cls,
        adata,
        metrics,
        mode,
        search_space,
        num_samples,
        max_epochs=max_epochs,
        scheduler=scheduler,
        searcher=searcher,
        seed=seed,
        resources=resources,
        name=experiment_name,
        logging_dir=logging_dir,
        scheduler_kwargs=scheduler_kwargs,
        searcher_kwargs=searcher_kwargs,
    )
    logging.info(f"Running autotune experiment {experiment.name}.")
    init(log_to_driver=False)
    experiment.result_grid = _configure_tuner(experiment).fit()
    return experiment
