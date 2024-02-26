from __future__ import annotations

from typing import Literal

from scvi._types import AnnOrMuData
from scvi.autotune._experiment import Experiment
from scvi.autotune._manager import TunerManager
from scvi.model.base import BaseModelClass


class ModelTuner:
    """``BETA`` Automated and scalable hyperparameter tuning for scvi-tools models.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.

    Parameters
    ----------
    model_cls
        A model class on which to tune hyperparameters. Must have a class property
        `_tunables` that defines tunable elements.

    Examples
    --------
    >>> import anndata
    >>> import scvi
    >>> adata = anndata.read_h5ad(path_to_h5ad)
    >>> model_cls = scvi.model.SCVI
    >>> model_cls.setup_anndata(adata)
    >>> tuner = scvi.autotune.ModelTuner(model_cls)
    >>> results = tuner.fit(adata, metric="validation_loss")

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/tuning/autotune_new_model`
    2. :doc:`/tutorials/notebooks/tuning/autotune_scvi`

    Lifecycle: beta
    """

    def __init__(self, model_cls: BaseModelClass):
        self.manager = TunerManager(model_cls)

    @property
    def manager(self) -> TunerManager:
        """The manager instance associated with the model class."""
        return self._manager

    @manager.setter
    def manager(self, value: TunerManager) -> None:
        if hasattr(self, "_manager"):
            raise AttributeError("Cannot reassign `manager`.")
        self._manager = value

    @property
    def history(self) -> dict[str, Experiment]:
        """A dictionary of completed tuning experiments."""
        return self.manager.history

    def info(self, extended_info: bool = False) -> None:
        self.manager.view_registry(extended_info=extended_info)

    def fit(
        self,
        adata: AnnOrMuData,
        metric: str,
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
        additional_metrics: list[str] | None = None,
        model_kwargs: dict | None = None,
        train_kwargs: dict | None = None,
        max_epochs: int | None = None,
        scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
        searcher: Literal["hyperopt", "random"] = "hyperopt",
        seed: int | None = None,
        resources: dict[str, float] | None = None,
        experiment_name: str | None = None,
        logging_dir: str | None = None,
        scheduler_kwargs: dict | None = None,
        searcher_kwargs: dict | None = None,
    ) -> None:
        """Run a specified hyperparameter sweep for the associated model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup with the
            associated model class.
        metric
            Primary metric to optimize.
        mode
            Optimization mode for the primary metric. One of ``"min"`` or ``"max"``.
        search_space
            Dictionary of hyperparameter names and their respective search spaces. See
            the `API <https://docs.ray.io/en/latest/tune/api/search_space.html>`_ for available
            search specifications.
        num_samples
            Total number of hyperparameter configurations to sample from the search space. Passed
            into :class:`~ray.tune.tune_config.TuneConfig`.
        additional_metrics
            Additional metrics to track during the experiment.
        model_kwargs
            Keyword arguments passed to the model class's constructor. Arguments must not overlap
            with those in ``search_space``.
        train_kwargs
            Keyword arguments passed to the model's ``train`` method. Arguments must not overlap
            with those in ``search_space``.

        max_epochs
            Maximum number of epochs to train each hyperparameter configuration.
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
            Random seed to use for the experiment.
        resources
            Dictionary of resources to allocate for each trial in the experiment. Available keys
            include:

            * ``"cpu"``: number of CPU cores
            * ``"gpu"``: number of GPUs
            * ``"memory"``: amount of memory

            If not provided, defaults to using all available resources. Note that fractional
            allocations are supported.
        experiment_name
            Name of the experiment, used for logging purposes. Defaults to a unique
            string formatted with the current timestamp and model class name.
        logging_dir
            Directory to store experiment logs. Defaults to a directory named ``"autotune"``
            in :attr:``scvi.settings.logging_dir``.
        scheduler_kwargs
            Additional keyword arguments to pass to the scheduler.
        searcher_kwargs
            Additional keyword arguments to pass to the search algorithm.
        """
        tuner, experiment = self.manager.configure_tuner(
            adata,
            metric,
            mode,
            search_space,
            num_samples,
            additional_metrics=additional_metrics,
            model_kwargs=model_kwargs,
            train_kwargs=train_kwargs,
            max_epochs=max_epochs,
            scheduler=scheduler,
            searcher=searcher,
            seed=seed,
            resources=resources,
            experiment_name=experiment_name,
            logging_dir=logging_dir,
            scheduler_kwargs=scheduler_kwargs,
            searcher_kwargs=searcher_kwargs,
        )
        return tuner.fit()

    def __repr__(self) -> str:
        return f"ModelTuner for {self.manager.model_cls.__name__}"
