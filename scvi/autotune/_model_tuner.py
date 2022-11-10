from typing import List, Optional

try:
    from ray import tune  # noqa
except ImportError:
    # ray tune dependency raises error in functions that use it
    pass

from scvi._compat import Literal
from scvi._types import AnnOrMuData
from scvi.autotune._manager import TunerManager
from scvi.model.base import BaseModelClass


class ModelTuner:
    """
    Automated and parallel hyperparameter searches with Ray Tune.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.
    Note: This API is in beta and is subject to change in future releases.

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
    >>> model_cls.setup_anndata(adata)
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
        num_samples: Optional[int] = None,
        use_defaults: bool = True,
        exclude: Optional[List[str]] = None,
        scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
        scheduler_kwargs: Optional[dict] = None,
        searcher: Literal["random", "grid", "hyperopt"] = "hyperopt",
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
    ) -> None:
        """
        Run a specified hyperparameter sweep over the model class.

        Does not require `adata` to be registered with the model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` object to use for
            training and validation.
        metric
            Primary metric to optimize over. If not provided, defaults to the model's
            validation loss.
        additional_metrics
            Additional metrics to report. If not provided, defaults to no metrics.
        search_space
            Dictionary of hyperparameter names and their respective search spaces.
            If not provided, defaults to the model's default search space. Options can
            be found in the model's documentation or with :func:`~scvi.autotune.ModelTuner.info`.
        num_samples
            Number of hyperparameter configurations to try. If not provided, defaults to
            an adaptive number of samples based on the available resources.
        use_defaults
            Whether to incorporate the model's default search space into the user-provided
            search space. Will default to `True` if `search_space` is not provided.
        exclude
            List of hyperparameter names in the model's default search space to exclude.
            Only used if `use_defaults` is `True`.
        scheduler
            Ray Tune scheduler to use. Supported options are:

            * ``"asha"``: Asynchronous Successive Halving Algorithm
            * ``"hyperband"``: Hyperband
            * ``"median"``: Median Stopping Rule
            * ``"pbt"``: Population Based Training
            * ``"fifo"``: First In First Out
        scheduler_kwargs
            Keyword arguments to pass to the scheduler.
        searcher
            Ray Tune searcher to use.
        searcher_kwargs
            Keyword arguments to pass to the searcher.
        reporter
            Whether to report progress using the Ray Tune reporter.
        resources
            Dictionary of maximum resources to allocate for the experiment. Defaults to
            using all available cpus and one gpu if available.
        setup_kwargs
            Keyword arguments to pass to the model's `setup_anndata` method.

        Returns
        -------
        None
        """
        tuner = self._manager.get_tuner(
            adata,
            metric=metric,
            additional_metrics=additional_metrics,
            search_space=search_space,
            num_samples=num_samples,
            use_defaults=use_defaults,
            exclude=exclude,
            scheduler=scheduler,
            scheduler_kwargs=scheduler_kwargs,
            searcher=searcher,
            searcher_kwargs=searcher_kwargs,
            reporter=reporter,
            resources=resources,
        )
        results = tuner.fit()
        return results

    def info(self, show_resources: bool = False):
        """View information about the current tuner."""
        self._manager._view_registry(show_resources=show_resources)
