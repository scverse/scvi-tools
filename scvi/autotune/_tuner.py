from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

from ._manager import TunerManager


class ModelTuner:
    """
    Automated and scalable hyperparameter tuning for scvi-tools models.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.
    Note: this API is in beta and is subject to change in future releases.

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
    """

    def __init__(self, model_cls: BaseModelClass):
        self._manager = TunerManager(model_cls)

    def fit(self, adata: AnnOrMuData, **kwargs) -> None:
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
            Additional metrics to track during the experiment. Defaults to `None`.
        search_space
            Dictionary of hyperparameter names and their respective search spaces
            provided as instantiated Ray Tune sample functions. Available
            hyperparameters can be viewed with :meth:`~scvi.autotune.ModelTuner.info`.
            Must be provided if `use_defaults` is `False`.
        use_defaults
            Whether to use the model class's default search space, which can be viewed
            with :meth:`~scvi.autotune.ModelTuner.info`. If `True` and `search_space` is
            provided, the two will be merged, giving priority to user-provided values.
            Defaults to `False`.
        num_samples
            Number of hyperparameter configurations to sample. Defaults to 10.
        max_epochs
            Maximum number of epochs to train each model configuration. Defaults to 100.
        scheduler
            Ray Tune scheduler to use. One of the following:

            * ``"asha"``: :class:`~ray.tune.schedulers.AsyncHyperBandScheduler` (default)
            * ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
            * ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
            * ``"pbt"``: :class:`~ray.tune.schedulers.PopulationBasedTraining`
            * ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`

            Note that that not all schedulers are compatible with all search algorithms.
            See Ray Tune `documentation <https://docs.ray.io/en/latest/tune/key-concepts.html#schedulers>`_
            for more details.
        scheduler_kwargs
            Keyword arguments to pass to the scheduler.
        searcher
            Ray Tune search algorithm to use. One of the following:

            * ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator` (default)
            * ``"hyperopt"``: :class:`~ray.tune.hyperopt.HyperOptSearch`
        searcher_kwargs
            Keyword arguments to pass to the search algorithm.
        reporter
            Whether to display progress with a Ray Tune reporter. Defaults to `True`.
            Depending on the execution environment, one of the following:

            * :class:`~ray.tune.CLIReporter` if running non-interactively
            * :class:`~ray.tune.JupyterNotebookReporter` if running interatively
        resources
            Dictionary of maximum resources to allocate for the experiment. Available
            keys include:

            * ``"cpu"``: number of CPU threads
            * ``"gpu"``: number of GPUs
            * ``"memory"``: amount of memory

            If not provided, defaults to using all available resources. Note that
            fractional allocations are supported.
        experiment_name
            Name of the experiment, used for logging purposes. Defaults to a unique
            string with the format `"tune_{model_cls}_{timestamp}"`.
        logging_dir
            Directory to store experiment logs. Defaults to a directory named `ray` in
            the current working directory.

        Returns
        -------
        :class:`~scvi.autotune.TuneAnalysis`
            A dataclass containing the results of the tuning experiment.
        """
        tuner, config = self._manager._get_tuner(adata, **kwargs)
        results = tuner.fit()
        return self._manager._get_analysis(results, config)

    def info(self, **kwargs) -> None:  # noqa: D102
        self._manager._view_registry(**kwargs)

    info.__doc__ = TunerManager._view_registry.__doc__
