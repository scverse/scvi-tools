from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

from ._manager import TunerManager


class ModelTuner:
    """
    Automated and parallel hyperparameter tuning with :ref:`~ray.tune`.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.
    Note: this API is in beta and is subject to change in future releases.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters.
        Currently supported model classes are:

        * :class:`~scvi.model.SCVI`

    Examples
    --------
    >>> import anndata
    >>> import scvi
    >>> adata = anndata.read_h5ad(path_to_h5ad)
    >>> model_cls = scvi.model.SCVI
    >>> model_cls.setup_anndata(adata)
    >>> tuner = scvi.autotune.ModelTuner(model_cls)
    >>> results = tuner.fit(adata, metric="validation_loss)
    """

    def __init__(self, model_cls: BaseModelClass):
        self._manager = TunerManager(model_cls)

    def fit(
        self,
        adata: AnnOrMuData,
        **kwargs,
    ) -> None:
        """
        Run a specified hyperparameter sweep for the associated model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup
            with the associated model class.
        metric
            The primary metric to optimize. If not provided, defaults to the model
            class's validation loss.
        additional_metrics
            Additional metrics to track during the experiment. If not provided, defaults
            to no other metrics.
        search_space
            Dictionary of hyperparameter names and their respective search spaces
            provided as instantiated Ray Tune sample functions. Available
            hyperparameters can be viewed with :meth:`~scvi.autotune.ModelTuner.info`.
            Must be provided if `use_defaults` is `False`.
        use_defaults
            Whether to use the model class's default search space, which can be viewed
            with :meth:`~scvi.autotune.ModelTuner.info`. If `True` and `search_space` is
            provided, the two will be merged, giving priority to user-provided values.
            Defaults to `True`.
        exclude
            List of hyperparameters to exclude from the default search space. If
            `use_defaults` is `False`, this argument is ignored.
        num_samples
            Number of hyperparameter configurations to sample.
        max_epochs
            Maximum number of epochs to train each model.
        scheduler
            Ray Tune scheduler to use. Supported options are:

            * ``"asha"``: :class:`~ray.tune.schedulers.ASHAScheduler`
            * ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
            * ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
            * ``"pbt"``: :class:`~ray.tune.schedulers.PopulationBasedTraining`
            * ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`
        scheduler_kwargs
            Keyword arguments to pass to the scheduler.
        searcher
            Ray Tune search algorithm to use. Supported options are:

            * ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`
            * ``"grid"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`
            * ``"hyperopt"``: :class:`~ray.tune.hyperopt.HyperOptSearch`
        searcher_kwargs
            Keyword arguments to pass to the search algorithm.
        reporter
            Whether to display progress with a Ray Tune reporter. Depending on the
            execution environment, will use one of the following reporters:

            * :class:`~ray.tune.CLIReporter` if running in a script
            * :class:`~ray.tune.JupyterNotebookReporter` if running in a notebook
        resources
            Dictionary of maximum resources to allocate for the experiment. Available
            keys include:

            * ``"cpu"``: maximum number of CPU threads to use
            * ``"gpu"``: maximum number of GPUs to use

            If not provided, defaults to using one CPU thread and one GPU if available.
        """
        tuner = self._manager._get_tuner(adata, **kwargs)
        results = tuner.fit()
        return results

    def info(self, show_resources: bool = False) -> None:
        """Display information about the associated model class."""
        self._manager._view_registry(show_resources=show_resources)
