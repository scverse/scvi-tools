from __future__ import annotations

from os.path import join
from typing import Any, Literal

from lightning.pytorch.callbacks import Callback
from lightning.pytorch.loggers import TensorBoardLogger
from ray.tune import ResultGrid, Tuner
from ray.tune.schedulers import TrialScheduler
from ray.tune.search import SearchAlgorithm

from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

_ASHA_DEFAULT_KWARGS = {
    "max_t": 100,
    "grace_period": 1,
    "reduction_factor": 2,
}


class AutotuneExperiment:
    """``BETA`` Track hyperparameter tuning experiments.

    Parameters
    ----------
    model_cls
        Model class on which to tune hyperparameters. Must implement a constructor and a ``train``
        method.
    adata
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup with
        ``model_cls``.
    metrics
        Either a single metric or a list of metrics to track during the experiment. If a list is
        provided, the primary metric will be the first element in the list.
    mode
        Optimization mode for the primary metric. One of ``"min"`` or ``"max"``.
    search_space
        Dictionary of hyperparameter names and their respective search spaces. See the
        `API <https://docs.ray.io/en/latest/tune/api/search_space.html>`_ for available search
        specifications. Must only contain the following top-level keys:

        * ``"model_args"``: parameters to pass to the model constructor.
        * ``"train_args"``: parameters to pass to the model's ``train`` method.

        Passed into :class:`~ray.tune.Tuner` as ``param_space``.
    num_samples
        Total number of hyperparameter configurations to sample from the search space. Passed into
        :class:`~ray.tune.tune_config.TuneConfig`.
    scheduler
        Ray Tune scheduler to use. One of the following:

        * ``"asha"``: :class:`~ray.tune.schedulers.AsyncHyperBandScheduler`
        * ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
        * ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
        * ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`

        Configured with reasonable defaults, which can be overridden with ``scheduler_kwargs``.
    searcher
        Ray Tune search algorithm to use. One of the following:

        * ``"hyperopt"``: :class:`~ray.tune.search.hyperopt.HyperOptSearch`
        * ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`

        Configured with reasonable defaults, which can be overridden with ``searcher_kwargs``.
    seed
        Random seed to use for the experiment. Propagated to :attr:`~scvi.settings.seed` and search
        algorithms. If not provided, defaults to :attr:`~scvi.settings.seed`.
    resources
        Dictionary of resources to allocate per trial in the experiment. Available keys include:

        * ``"cpu"``: number of CPU cores
        * ``"gpu"``: number of GPUs
        * ``"memory"``: amount of memory

        Passed into :func:`~ray.tune.with_resources`.
    name
        Name of the experiment, used for logging purposes. Defaults to a unique ID.
    logging_dir
        Base directory to store experiment logs. Defaults to :attr:``scvi.settings.logging_dir``.
    scheduler_kwargs
        Additional keyword arguments to pass to the scheduler.
    searcher_kwargs
        Additional keyword arguments to pass to the search algorithm.

    Notes
    -----
    Lifecycle: beta

    See Also
    --------
    :func:`~scvi.autotune.run_autotune`
    """

    def __init__(
        self,
        model_cls: BaseModelClass,
        adata: AnnOrMuData,
        metrics: str | list[str],
        mode: Literal["min", "max"],
        search_space: dict[str, dict[Literal["model_args", "train_args"], dict[str, Any]]],
        num_samples: int,
        scheduler: Literal["asha", "hyperband", "median", "fifo"] = "asha",
        searcher: Literal["hyperopt", "random"] = "hyperopt",
        seed: int | None = None,
        resources: dict[Literal["cpu", "gpu", "memory"], float] | None = None,
        name: str | None = None,
        logging_dir: str | None = None,
        scheduler_kwargs: dict | None = None,
        searcher_kwargs: dict | None = None,
    ) -> None:
        self.model_cls = model_cls
        self.adata = adata
        self.metrics = metrics
        self.mode = mode
        self.search_space = search_space
        self.num_samples = num_samples
        self.seed = seed
        self.scheduler_kwargs = scheduler_kwargs
        self.searcher_kwargs = searcher_kwargs
        self.scheduler = scheduler
        self.searcher = searcher
        self.resources = resources
        self.name = name
        self.logging_dir = logging_dir

    @property
    def id(self) -> str:
        """Unique identifier for the experiment."""
        from uuid import uuid4

        if not hasattr(self, "_id"):
            self._id = str(uuid4())
        return self._id

    @property
    def model_cls(self) -> BaseModelClass:
        """Model class on which to tune hyperparameters."""
        return self._model_cls

    @model_cls.setter
    def model_cls(self, value: BaseModelClass) -> None:
        if hasattr(self, "_model_cls"):
            raise AttributeError("Cannot reassign `model_cls`")
        self._model_cls = value

    @property
    def adata(self) -> AnnOrMuData:
        """:class:`~anndata.AnnData` or :class:`~mudata.MuData` for the experiment."""
        return self._adata

    @adata.setter
    def adata(self, value: AnnOrMuData) -> None:
        from scvi.data._constants import _SETUP_ARGS_KEY, _SETUP_METHOD_NAME

        if hasattr(self, "_adata"):
            raise AttributeError("Cannot reassign `adata`")

        data_manager = self.model_cls._get_most_recent_anndata_manager(value, required=True)
        self._adata = value
        self._setup_method_name = data_manager._registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        self._setup_method_args = data_manager._get_setup_method_args().get(_SETUP_ARGS_KEY, {})

    @property
    def setup_method_name(self) -> str:
        """Either ``"setup_anndata"`` or ``"setup_mudata"``."""
        if not hasattr(self, "_setup_method_name"):
            raise AttributeError("`setup_method_name` not yet available.")
        return self._setup_method_name

    @property
    def setup_method_args(self) -> dict[str, Any]:
        """Keyword arguments for the setup method."""
        if not hasattr(self, "_setup_method_args"):
            raise AttributeError("`setup_method_args` not yet available.")
        return self._setup_method_args

    @property
    def metrics(self) -> list[str]:
        """Metrics to track during the experiment."""
        if not hasattr(self, "_metrics"):
            raise AttributeError("`metrics` not yet available.")
        return self._metrics

    @metrics.setter
    def metrics(self, value: str | list[str]) -> None:
        if hasattr(self, "_metrics"):
            raise AttributeError("Cannot reassign `metrics`")
        elif value is None:
            raise ValueError("`metrics` must not be `None`")
        elif not isinstance(value, str) and not isinstance(value, list):
            raise TypeError("`metrics` must be a string or a list of strings")
        elif isinstance(value, str):
            value = [value]
        elif isinstance(value, list) and len(value) < 1:
            raise ValueError("`metrics` must not be empty")
        self._metrics = value

    @property
    def mode(self) -> Literal["min", "max"]:
        """Optimization mode for the primary metric."""
        if not hasattr(self, "_mode"):
            raise AttributeError("`mode` not yet available.")
        return self._mode

    @mode.setter
    def mode(self, value: Literal["min", "max"]) -> None:
        if hasattr(self, "_mode"):
            raise AttributeError("Cannot reassign `mode`")
        elif value not in ["min", "max"]:
            raise ValueError("`mode` must be either 'min' or 'max'")
        self._mode = value

    @property
    def search_space(self) -> dict[str, dict[Literal["model_args", "train_args"], dict[str, Any]]]:
        """Search space for hyperparameters."""
        if not hasattr(self, "_search_space"):
            raise AttributeError("`search_space` not yet available.")
        return self._search_space

    @search_space.setter
    def search_space(self, value: dict[str, Any]) -> None:
        if hasattr(self, "_search_space"):
            raise AttributeError("Cannot reassign `search_space`")
        elif not isinstance(value, dict):
            raise TypeError("`search_space` must be a dictionary")
        elif len(value) == 0:
            raise ValueError("`search_space` must not be empty")
        elif any(key not in ["model_args", "train_args"] for key in value.keys()):
            raise KeyError("`search_space` must only contain 'model_args' and 'train_args' keys")

        self._search_space = value

    @property
    def num_samples(self) -> int:
        """Total number of hyperparameter configurations to sample."""
        if not hasattr(self, "_num_samples"):
            raise AttributeError("`num_samples` not yet available.")
        return self._num_samples

    @num_samples.setter
    def num_samples(self, value: int | None) -> None:
        if hasattr(self, "_num_samples"):
            raise AttributeError("Cannot reassign `num_samples`")
        elif not isinstance(value, int):
            raise TypeError("`num_samples` must be an integer")
        self._num_samples = value

    @property
    def scheduler(self) -> TrialScheduler:
        """Ray Tune scheduler to use."""
        if not hasattr(self, "_scheduler"):
            raise AttributeError("`scheduler` not yet available.")
        return self._scheduler

    @scheduler.setter
    def scheduler(self, value: Literal["asha", "hyperband", "median", "fifo"]) -> None:
        from ray.tune.schedulers import (
            AsyncHyperBandScheduler,
            FIFOScheduler,
            HyperBandScheduler,
            MedianStoppingRule,
        )

        if hasattr(self, "_scheduler"):
            raise AttributeError("Cannot reassign `scheduler`")
        elif not isinstance(value, str):
            raise TypeError("`scheduler` must be a string")
        elif value not in ["asha", "hyperband", "median", "fifo"]:
            raise ValueError("`scheduler` must be one of 'asha', 'hyperband', 'median', 'fifo'")

        kwargs = {
            "metric": self.metrics[0],
            "mode": self.mode,
        }
        if value == "asha":
            kwargs.update(_ASHA_DEFAULT_KWARGS)
            scheduler_cls = AsyncHyperBandScheduler
        elif value == "hyperband":
            scheduler_cls = HyperBandScheduler
        elif value == "median":
            scheduler_cls = MedianStoppingRule
        else:
            kwargs = {}
            scheduler_cls = FIFOScheduler

        kwargs.update(self.scheduler_kwargs)
        self._scheduler = scheduler_cls(**kwargs)

    @property
    def scheduler_kwargs(self) -> dict[str, Any]:
        """Additional keyword arguments to pass to the scheduler."""
        if not hasattr(self, "_scheduler_kwargs"):
            raise AttributeError("`scheduler_kwargs` not yet available.")
        return self._scheduler_kwargs

    @scheduler_kwargs.setter
    def scheduler_kwargs(self, value: dict | None) -> None:
        if hasattr(self, "_scheduler_kwargs"):
            raise AttributeError("Cannot reassign `scheduler_kwargs`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`scheduler_kwargs` must be a dictionary")
        self._scheduler_kwargs = value or {}

    @property
    def searcher(self) -> SearchAlgorithm:
        """Ray Tune search algorithm to use."""
        if not hasattr(self, "_searcher"):
            raise AttributeError("`searcher` not yet available.")
        return self._searcher

    @searcher.setter
    def searcher(self, value: Literal["hyperopt", "random"]) -> None:
        from ray.tune.search import BasicVariantGenerator
        from ray.tune.search.hyperopt import HyperOptSearch

        if hasattr(self, "_searcher"):
            raise AttributeError("Cannot reassign `searcher`")
        elif not isinstance(value, str):
            raise TypeError("`searcher` must be a string")
        elif value not in ["hyperopt", "random"]:
            raise ValueError("`searcher` must be one of 'hyperopt', 'random'")

        if value == "random":
            kwargs = {"random_state": self.seed}
            searcher_cls = BasicVariantGenerator
        else:
            kwargs = {
                "metric": self.metrics[0],
                "mode": self.mode,
                "random_state_seed": self.seed,
            }
            searcher_cls = HyperOptSearch

        kwargs.update(self.searcher_kwargs)
        self._searcher = searcher_cls(**kwargs)

    @property
    def searcher_kwargs(self) -> dict[str, Any]:
        """Additional keyword arguments to pass to the search algorithm."""
        if not hasattr(self, "_searcher_kwargs"):
            raise AttributeError("`searcher_kwargs` not yet available.")
        return self._searcher_kwargs

    @searcher_kwargs.setter
    def searcher_kwargs(self, value: dict | None) -> None:
        if hasattr(self, "_searcher_kwargs"):
            raise AttributeError("Cannot reassign `searcher_kwargs`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`searcher_kwargs` must be a dictionary")
        self._searcher_kwargs = value or {}

    @property
    def seed(self) -> int | None:
        """Random seed to use for the experiment."""
        if not hasattr(self, "_seed"):
            raise AttributeError("`seed` not yet available.")
        return self._seed

    @seed.setter
    def seed(self, value: int | None) -> None:
        from scvi import settings

        if hasattr(self, "_seed"):
            raise AttributeError("Cannot reassign `seed`")
        elif value is not None and not isinstance(value, int):
            raise TypeError("`seed` must be an integer")
        self._seed = value or settings.seed

    @property
    def resources(self) -> dict[Literal["cpu", "gpu", "memory"], float] | None:
        """Resources to allocate per trial in the experiment."""
        if not hasattr(self, "_resources"):
            raise AttributeError("`resources` not yet available.")
        return self._resources

    @resources.setter
    def resources(self, value: dict[Literal["cpu", "gpu", "memory"], float] | None) -> None:
        if hasattr(self, "_resources"):
            raise AttributeError("Cannot reassign `resources`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`resources` must be a dictionary or `None`")
        self._resources = value or {}

    @property
    def name(self) -> str:
        """Name of the experiment."""
        if not hasattr(self, "_name"):
            raise AttributeError("`name` not yet available.")
        return self._name

    @name.setter
    def name(self, value: str | None) -> None:
        if hasattr(self, "_name"):
            raise AttributeError("Cannot reassign `name`")
        elif value is not None and not isinstance(value, str):
            raise TypeError("`name` must be a string or `None`")

        if value is None:
            default = f"{self._model_cls.__name__.lower()}_"
            default += self.id
        self._name = value or default

    @property
    def logging_dir(self) -> str:
        """Base directory to store experiment logs."""
        if not hasattr(self, "_logging_dir"):
            raise AttributeError("`logging_dir` not yet available.")
        return self._logging_dir

    @logging_dir.setter
    def logging_dir(self, value: str | None) -> None:
        from scvi import settings

        if hasattr(self, "_logging_dir"):
            raise AttributeError("Cannot reassign `logging_dir`")
        elif value is not None and not isinstance(value, str):
            raise TypeError("`logging_dir` must be a string")
        self._logging_dir = value or join(settings.logging_dir, self.name)

    @property
    def metrics_callback(self) -> Callback:
        from ray.tune.integration.pytorch_lightning import TuneReportCheckpointCallback

        if not hasattr(self, "_metrics"):
            raise AttributeError("`metrics_callback` not yet available.")

        # this is to get around lightning import changes
        callback_cls = type(
            "_TuneReportCheckpointCallback",
            (TuneReportCheckpointCallback, Callback),
            {},
        )
        return callback_cls(metrics=self.metrics, on="validation_end", save_checkpoints=False)

    @property
    def result_grid(self) -> ResultGrid:
        """Result grid for the experiment."""
        if not hasattr(self, "_result_grid"):
            raise AttributeError("`result_grid` not yet available.")
        return self._result_grid

    @result_grid.setter
    def result_grid(self, value: ResultGrid) -> None:
        if hasattr(self, "_result_grid"):
            raise AttributeError("Cannot reassign `result_grid`")
        self._result_grid = value

    def __repr__(self) -> str:
        return f"Experiment {self.name}"

    def get_tuner(self) -> Tuner:
        """Configure a :class:`~ray.tune.Tuner` from this experiment."""
        from ray.train import RunConfig
        from ray.tune import with_parameters, with_resources
        from ray.tune.tune_config import TuneConfig

        trainable = with_parameters(_trainable, experiment=self)
        trainable = with_resources(trainable, resources=self.resources)

        tune_config = TuneConfig(
            scheduler=self.scheduler,
            search_alg=self.searcher,
            num_samples=self.num_samples,
        )
        run_config = RunConfig(
            name=self.name,
            storage_path=self.logging_dir,
            local_dir=self.logging_dir,
            log_to_file=True,
            verbose=1,
        )
        return Tuner(
            trainable=trainable,
            param_space=self.search_space,
            tune_config=tune_config,
            run_config=run_config,
        )

    def get_logger(self, trial_name: str) -> TensorBoardLogger:
        """Configure TensorBoard logger for a trial in this experiment."""
        return TensorBoardLogger(join(self.logging_dir, f"{trial_name}_tensorboard"))


def _trainable(
    param_sample: dict[str, dict[Literal["model_args", "train_args"], dict[str, Any]]],
    experiment: AutotuneExperiment,
) -> None:
    """Implements a Ray Tune trainable function for an :class:`~scvi.autotune.AutotuneExperiment`.

    Setup on the :class:`~anndata.AnnData` or :class:`~mudata.MuData` has to be performed since Ray
    opens a new process per trial and thus the initial setup on the main process is not
    transferred.

    Parameters
    ----------
    param_sample
        Hyperparameter configuration on which to train. Note: this is different from
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
    from ray.train import get_context

    from scvi import settings

    model_args, train_args = param_sample.get("model_args", {}), param_sample.get("train_args", {})
    train_args = {
        "accelerator": "auto",
        "devices": "auto",
        "check_val_every_n_epoch": 1,
        "enable_progress_bar": False,
        "logger": experiment.get_logger(get_context().get_trial_name()),
        "callbacks": [experiment.metrics_callback],
        **train_args,
    }

    settings.seed = experiment.seed
    getattr(experiment.model_cls, experiment.setup_method_name)(
        experiment.adata,
        **experiment.setup_method_args,
    )
    model = experiment.model_cls(experiment.adata, **model_args)
    model.train(**train_args)
