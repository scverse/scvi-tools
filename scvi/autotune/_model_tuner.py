from typing import List, Literal, Optional

from scvi._types import AnnOrMuData
from scvi.autotune._manager import TunerManager
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

    Examples
    --------
    >>> import scvi
    >>> model_cls = scvi.model.SCVI
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> tuner = ModelTuner(model_cls)
    >>> results = tuner.fit(adata, metric="validation_loss")
    """

    def __init__(
        self,
        model_cls: BaseModelClass,
    ):
        self._manager = TunerManager(model_cls)

    def fit(
        self,
        adata: AnnOrMuData,
        *,
        metric: Optional[str] = None,
        additional_metrics: Optional[List[str]] = None,
        search_space: Optional[dict] = None,
        use_defaults: bool = True,
        exclude: Optional[dict] = None,
        scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
        scheduler_kwargs: Optional[dict] = None,
        searcher: Literal["random", "grid", "hyperopt"] = "hyperopt",
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
        setup_kwargs: Optional[dict] = None,
    ) -> None:
        """
        Run a specified hyperparameter sweep over the model class with the given
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` object and search
        space.

        Does not require `adata` to be registered with the model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` object to use for
            training and validation.

        Returns
        -------

        """
        tuner = self._manager.get_tuner(
            adata,
            metric=metric,
            additional_metrics=additional_metrics,
            search_space=search_space,
            use_defaults=use_defaults,
            exclude=exclude,
            scheduler=scheduler,
            scheduler_kwargs=scheduler_kwargs,
            searcher=searcher,
            searcher_kwargs=searcher_kwargs,
            reporter=reporter,
            resources=resources,
            setup_kwargs=setup_kwargs,
        )
        results = tuner.fit()
        return results

    def info(self):
        self._manager.view_registry()
