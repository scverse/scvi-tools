from typing import Any, Literal

from ray import tune
from ray.tune import ResultGrid
from ray.tune.schedulers import TrialScheduler
from ray.tune.search import SearchAlgorithm

from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

_ASHA_DEFAULT_KWARGS = {
    "max_t": 100,
    "grace_period": 1,
    "reduction_factor": 2,
}


class Experiment:
    def __init__(
        self,
        model_cls: BaseModelClass,
        adata: AnnOrMuData,
        metrics: str | list[str],
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
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
        self.model_cls = model_cls
        self.adata = adata

        self.metrics = metrics
        self.mode = mode
        self.search_space = search_space
        self.num_samples = num_samples

        self.model_kwargs = model_kwargs
        self.train_kwargs = train_kwargs
        self.max_epochs = max_epochs

        self.scheduler_kwargs = scheduler_kwargs
        self.searcher_kwargs = searcher_kwargs
        self.scheduler = scheduler
        self.searcher = searcher
        self.validate_scheduler_searcher_compat(scheduler, searcher)

        self.seed = seed
        self.resources = resources
        self.experiment_name = experiment_name
        self.logging_dir = logging_dir

    @property
    def id(self) -> str:
        if not hasattr(self, "_id"):
            from uuid import uuid4

            self._id = str(uuid4())
        return self._id

    @property
    def adata(self) -> AnnOrMuData:
        return self._adata

    @adata.setter
    def adata(self, value: AnnOrMuData) -> None:
        from scvi.data._constants import _SETUP_ARGS_KEY, _SETUP_METHOD_NAME

        if hasattr(self, "_adata"):
            raise AttributeError("Cannot reassign `adata`")

        data_manager = self.model_cls._get_most_recent_anndata_manager(value, required=True)
        self._adata = value
        self._setup_method_name = data_manager._registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        self._setup_args = data_manager._get_setup_method_args().get(_SETUP_ARGS_KEY, {})

    @property
    def setup_method_name(self) -> str:
        return self._setup_method_name

    @property
    def setup_args(self) -> dict[str, Any]:
        return self._setup_args

    @property
    def model_cls(self) -> BaseModelClass:
        return self._model_cls

    @model_cls.setter
    def model_cls(self, value: BaseModelClass) -> None:
        if hasattr(self, "_model_cls"):
            raise AttributeError("Cannot reassign `model_cls`")
        self._model_cls = value

    @property
    def metrics(self) -> list[str]:
        return self._metrics

    @metrics.setter
    def metrics(self, value: str | list[str]) -> None:
        if hasattr(self, "_metrics"):
            raise AttributeError("Cannot reassign `metrics`")
        elif not isinstance(value, str) and not isinstance(value, list):
            raise TypeError("`metrics` must be a string or a list of strings")
        elif isinstance(value, str):
            value = [value]
        self._metrics = value

    @property
    def mode(self) -> Literal["min", "max"]:
        return self._mode

    @mode.setter
    def mode(self, value: Literal["min", "max"]) -> None:
        if hasattr(self, "_mode"):
            raise AttributeError("Cannot reassign `mode`")
        elif value not in ["min", "max"]:
            raise ValueError("`mode` must be either 'min' or 'max'")
        self._mode = value

    @property
    def search_space(self) -> dict[str, callable]:
        return self._search_space

    @search_space.setter
    def search_space(self, value: dict[str, callable]) -> None:
        if hasattr(self, "_search_space"):
            raise AttributeError("Cannot reassign `search_space`")
        elif not isinstance(value, dict):
            raise TypeError("`search_space` must be a dictionary")
        elif len(value) == 0:
            raise ValueError("`search_space` must not be empty")

        self._search_space = value

    @property
    def num_samples(self) -> int:
        return self._num_samples

    @num_samples.setter
    def num_samples(self, value: int | None) -> None:
        if hasattr(self, "_num_samples"):
            raise AttributeError("Cannot reassign `num_samples`")
        elif not isinstance(value, int):
            raise TypeError("`num_samples` must be an integer")
        self._num_samples = value

    @property
    def model_kwargs(self) -> dict[str, Any]:
        return self._model_kwargs

    @model_kwargs.setter
    def model_kwargs(self, value: dict | None) -> None:
        if hasattr(self, "_model_kwargs"):
            raise AttributeError("Cannot reassign `model_kwargs`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`model_kwargs` must be a dictionary")

        value = value or {}
        for param in value:
            if param in self.search_space:
                raise ValueError(
                    f"`model_kwargs` contains a parameter `{param}` that is also in `search_space`."
                )
        self._model_kwargs = value

    @property
    def train_kwargs(self) -> dict[str, Any]:
        return self._train_kwargs

    @train_kwargs.setter
    def train_kwargs(self, value: dict | None) -> None:
        if hasattr(self, "_train_kwargs"):
            raise AttributeError("Cannot reassign `train_kwargs`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`train_kwargs` must be a dictionary")

        value = value or {}
        for param in value:
            if param in self.search_space:
                raise ValueError(
                    f"`train_kwargs` contains a parameter `{param}` that is also in `search_space`."
                )
        self._train_kwargs = value

    @property
    def max_epochs(self) -> int | None:
        return self._max_epochs

    @max_epochs.setter
    def max_epochs(self, value: int | None) -> None:
        if hasattr(self, "_max_epochs"):
            raise AttributeError("Cannot reassign `max_epochs`")
        elif value is not None and not isinstance(value, int):
            raise TypeError("`max_epochs` must be an integer")
        self._max_epochs = value

    @property
    def scheduler(self) -> TrialScheduler:
        return self._scheduler

    @scheduler.setter
    def scheduler(self, value: Literal["asha", "hyperband", "median", "pbt", "fifo"]) -> None:
        if hasattr(self, "_scheduler"):
            raise AttributeError("Cannot reassign `scheduler`")
        elif not isinstance(value, str):
            raise TypeError("`scheduler` must be a string")
        elif value not in ["asha", "hyperband", "median", "pbt", "fifo"]:
            raise ValueError(
                "`scheduler` must be one of 'asha', 'hyperband', 'median', 'pbt', 'fifo'"
            )

        kwargs = {
            "metric": self.metrics[0],
            "mode": self.mode,
        }
        if value == "asha":
            from ray.tune.schedulers import AsyncHyperBandScheduler

            kwargs.update(_ASHA_DEFAULT_KWARGS)
            scheduler_cls = AsyncHyperBandScheduler
        elif value == "hyperband":
            from ray.tune.schedulers import HyperBandScheduler

            scheduler_cls = HyperBandScheduler
        elif value == "median":
            from ray.tune.schedulers import MedianStoppingRule

            scheduler_cls = MedianStoppingRule
        elif value == "pbt":
            from ray.tune.schedulers import PopulationBasedTraining

            scheduler_cls = PopulationBasedTraining
        else:
            from ray.tune.schedulers import FIFOScheduler

            kwargs = {}
            scheduler_cls = FIFOScheduler

        kwargs.update(self.scheduler_kwargs)
        self._scheduler = scheduler_cls(**kwargs)

    @property
    def scheduler_kwargs(self) -> dict[str, Any]:
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
        return self._searcher

    @searcher.setter
    def searcher(self, value: Literal["hyperopt", "random"]) -> None:
        if hasattr(self, "_searcher"):
            raise AttributeError("Cannot reassign `searcher`")
        elif not isinstance(value, str):
            raise TypeError("`searcher` must be a string")
        elif value not in ["hyperopt", "random"]:
            raise ValueError("`searcher` must be one of 'hyperopt', 'random'")

        kwargs = {
            "metric": self.metrics[0],
            "mode": self.mode,
        }
        if value == "random":
            from ray.tune.search import BasicVariantGenerator

            kwargs = {}
            searcher_cls = BasicVariantGenerator
        else:
            # tune does not import hyperopt by default
            tune.search.SEARCH_ALG_IMPORT[value]()
            searcher_cls = tune.search.hyperopt.HyperOptSearch

        kwargs.update(self.searcher_kwargs)
        self._searcher = searcher_cls(**kwargs)

    @property
    def searcher_kwargs(self) -> dict[str, Any]:
        return self._searcher_kwargs

    @searcher_kwargs.setter
    def searcher_kwargs(self, value: dict | None) -> None:
        if hasattr(self, "_searcher_kwargs"):
            raise AttributeError("Cannot reassign `searcher_kwargs`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`searcher_kwargs` must be a dictionary")
        self._searcher_kwargs = value or {}

    def validate_scheduler_searcher_compat(self, scheduler: str, searcher: str) -> None:
        if scheduler not in ["asha", "median", "hyperband"] and searcher not in ["random"]:
            raise ValueError(
                f"`scheduler={scheduler}` is incompatible with `searcher={searcher}`. Please see "
                f"https://docs.ray.io/en/latest/tune/key-concepts.html#tune-schedulers for more "
                "information."
            )

    @property
    def seed(self) -> int | None:
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
    def resources(self) -> dict[str, float]:
        return self._resources

    @resources.setter
    def resources(self, value: dict[str, float] | None) -> None:
        from math import ceil

        if hasattr(self, "_resources"):
            raise AttributeError("Cannot reassign `resources`")
        elif value is not None and not isinstance(value, dict):
            raise TypeError("`resources` must be a dictionary")
        self._resources = value or {}

        self._accelerator = "gpu" if self.resources.get("gpu", 0) > 0 else "cpu"
        self._devices = self.resources.get(self.accelerator, 1)
        if isinstance(self._devices, float):
            self._devices = ceil(self._devices)

    @property
    def accelerator(self) -> str:
        return self._accelerator

    @property
    def devices(self) -> int:
        return self._devices

    @property
    def experiment_name(self) -> str:
        return self._experiment_name

    @experiment_name.setter
    def experiment_name(self, value: str | None) -> None:
        if hasattr(self, "_experiment_name"):
            raise AttributeError("Cannot reassign `experiment_name`")
        elif value is not None and not isinstance(value, str):
            raise TypeError("`experiment_name` must be a string")

        if value is None:
            default = f"_{self._model_cls.__name__.lower()}"
            default += self.id
        self._experiment_name = value or default

    @property
    def logging_dir(self) -> str:
        return self._logging_dir

    @logging_dir.setter
    def logging_dir(self, value: str | None) -> None:
        from scvi import settings

        if hasattr(self, "_logging_dir"):
            raise AttributeError("Cannot reassign `logging_dir`")
        elif value is not None and not isinstance(value, str):
            raise TypeError("`logging_dir` must be a string")
        self._logging_dir = value or settings.logging_dir

    @property
    def result_grid(self) -> ResultGrid:
        if not hasattr(self, "_result_grid"):
            raise AttributeError("`result_grid` not yet available.")
        return self._result_grid

    @result_grid.setter
    def result_grid(self, value: ResultGrid) -> None:
        if hasattr(self, "_result_grid"):
            raise AttributeError("Cannot reassign `result_grid`")
        self._result_grid = value

    def __repr__(self) -> str:
        return f"Experiment {self.experiment_name} for {self.model_cls.__name__}"
