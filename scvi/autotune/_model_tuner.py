from typing import Any, Optional

from scvi._types import AnnOrMuData
from scvi.autotune import TunerManager
from scvi.model.base import BaseModelClass


class ModelTuner:
    """
    Automated and parallel hyperparameter searches with Ray Tune.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters.
        Currently supports one of the following:

        * :class:`~scvi.model.SCVI`
    adata
        AnnData/MuData object that has been registered via the model's ``setup_anndata``
        or ``setup_mudata`` method.
    """

    def __init__(
        self,
        model_cls: BaseModelClass,
    ):
        self._manager = TunerManager(model_cls)  # safe to import tune after this

    def fit(
        self,
        adata: AnnOrMuData,
        search_config: Optional[dict] = None,
        use_defaults: bool = True,
        exclude: Optional[dict] = None,
    ) -> None:
        from ray import tune

        search_config = self._manager.validate_search_space(
            search_config,
            use_defaults=use_defaults,
            exclude=exclude,
        )
        reporter = tune.CLIReporter(
            parameter_columns=None,
            metric_columns=None,
        )
        print(reporter)

    @staticmethod
    def _fit(
        model_cls: BaseModelClass,
        adata: AnnOrMuData,
        search_config: dict,
        scheduler: Any,
        reporter: Any,
    ) -> None:
        pass
        # from ray import tune

        # tuner = tune.Tuner(
        #     tune.Trainable,
        #     tune_config=None,
        #     run_config=None,
        #     param_space=None,
        # )
        # results = tuner.fit()
