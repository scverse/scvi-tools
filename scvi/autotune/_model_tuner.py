import logging
from typing import Any, Literal, Optional

from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

from ._manager import TunerManager

logger = logging.getLogger(__name__)


class ModelTuner:
    """
    Automated and parallel hyperparameter searches with Ray Tune. Wraps a
    :class:`~ray.tune.Tuner` instance that is attached to a scvi-tools model
    class.

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
        adata: AnnOrMuData,
    ):
        self._manager = TunerManager(model_cls, adata)

    def fit(
        self,
        metric: str,
        search_config: Optional[dict] = None,
        use_defaults: bool = True,
        scheduler: Literal["asha"] = "asha",
        resources_per_trial: Optional[dict] = None,
    ) -> None:
        # from ray import tune

        search_config = self._manager.validate_search_space(
            search_config, use_defaults=use_defaults
        )
        # scheduler = None
        # reporter = None
        # search = None

        return self._fit()

    @staticmethod
    def _fit(
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
