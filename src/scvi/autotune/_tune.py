from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi.autotune._experiment import AutotuneExperiment

if TYPE_CHECKING:
    from typing import Any, Literal

    from lightning.pytorch import LightningDataModule

    from scvi._types import AnnOrMuData
    from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


def run_autotune(
    model_cls: BaseModelClass,
    data: AnnOrMuData | LightningDataModule,
    metrics: str | list[str],
    mode: Literal["min", "max"],
    search_space: dict[str, dict[Literal["model_params", "train_params"], dict[str, Any]]],
    num_samples: int,
    scheduler: Literal["asha", "hyperband", "median", "fifo"] = "asha",
    searcher: Literal["hyperopt", "random"] = "hyperopt",
    seed: int | None = None,
    resources: dict[Literal["cpu", "gpu", "memory"], float] | None = None,
    experiment_name: str | None = None,
    logging_dir: str | None = None,
    scheduler_kwargs: dict | None = None,
    searcher_kwargs: dict | None = None,
    scib_stage: str | None = "train_end",
    scib_subsample_rows: int | None = 5000,
    scib_indices_list: list | None = None,
    log_to_driver: bool = False,
    local_mode: bool = False,
    ignore_reinit_error: bool = False,
) -> AutotuneExperiment:
    """``BETA`` Run a hyperparameter sweep.

    Parameters
    ----------
    model_cls
        Model class on which to tune hyperparameters.
    data
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup with
        ``model_cls`` or a :class:`~lightning.pytorch.core.LightningDataModule` (``EXPERIMENTAL``).
    metrics
        Either a single metric or a list of metrics to track during the experiment. If a list is
        provided, the primary metric will be the first element in the list.
    mode
        Optimization mode for the primary metric. One of ``"min"`` or ``"max"``.
    search_space
        Dictionary of hyperparameter names and their respective search spaces. See the
        `API <https://docs.ray.io/en/latest/tune/api/search_space.html>`_ for available search
        specifications. Must only contain the following top-level keys:

        * ``"model_params"``: parameters to pass to the model constructor.
        * ``"train_params"``: parameters to pass to the model's ``train`` method.

        Passed into :class:`~ray.tune.Tuner` as ``param_space``.
    num_samples
        Total number of hyperparameter configurations to sample from the search space. Passed into
        :class:`~ray.tune.tune_config.TuneConfig`.
    scheduler
        Ray Tune scheduler to use. One of the following:

        - ``"asha"``: :class:`~ray.tune.schedulers.AsyncHyperBandScheduler`
        - ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
        - ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
        - ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`

        Configured with reasonable defaults, which can be overridden with ``scheduler_kwargs``.
    searcher
        Ray Tune search algorithm to use. One of the following:

        - ``"hyperopt"``: :class:`~ray.tune.search.hyperopt.HyperOptSearch`
        - ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`

        Configured with reasonable defaults, which can be overridden with ``searcher_kwargs``.
    seed
        Random seed to use for the experiment. Propagated to :attr:`~scvi.settings.seed` and
        search algorithms. If not provided, defaults to :attr:`~scvi.settings.seed`.
    resources
        Dictionary of resources to allocate per trial in the experiment. Available keys
        include:

        - ``"cpu"``: number of CPU cores
        - ``"gpu"``: number of GPUs
        - ``"memory"``: amount of memory

        Passed into :func:`~ray.tune.with_resources`.
    experiment_name
        Name of the experiment, used for logging purposes. Defaults to a unique ID concatenated
        to the model class name.
    logging_dir
        Base directory to store experiment logs. Defaults to :attr:`~scvi.settings.logging_dir`.
    scheduler_kwargs
        Additional keyword arguments to pass to the scheduler.
    searcher_kwargs
        Additional keyword arguments to pass to the search algorithm.
    scib_stage
        Used when performing scib-metrics tune, select whether to perform on "validation_end"
        or "train_end" (default).
    scib_subsample_rows
        Used when performing scib-metrics tune, select number of rows to subsample (5000 default).
        This is important to save computation time
    scib_indices_list
        If not empty will be used to select the indices to calc the scib metric on, otherwise will
        use the random indices selection in size of scib_subsample_rows
    ignore_reinit_error
        If true, Ray suppresses errors from calling
        ray.init() a second time. Ray won't be restarted.
    local_mode
        Deprecated: consider using the Ray Debugger instead.
    log_to_driver
        If true, the output from all of the worker
        processes on all nodes will be directed to the driver.

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
        model_cls=model_cls,
        data=data,
        metrics=metrics,
        mode=mode,
        search_space=search_space,
        num_samples=num_samples,
        scheduler=scheduler,
        searcher=searcher,
        seed=seed,
        resources=resources,
        name=experiment_name,
        logging_dir=logging_dir,
        scheduler_kwargs=scheduler_kwargs,
        searcher_kwargs=searcher_kwargs,
        scib_stage=scib_stage,
        scib_subsample_rows=scib_subsample_rows,
        scib_indices_list=scib_indices_list,
    )
    logger.info(f"Running autotune experiment {experiment.name}.")
    init(
        log_to_driver=log_to_driver, ignore_reinit_error=ignore_reinit_error, local_mode=local_mode
    )
    experiment.result_grid = experiment.get_tuner().fit()
    return experiment
