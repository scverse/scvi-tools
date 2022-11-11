from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass


class ModelTuner:
    """
    Automated and parallel hyperparameter tuning with Ray Tune.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model
    class. Note: this API is in beta and is subject to change in future
    releases.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune
        hyperparameters. Currently supported model classes are:

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
        pass

    def fit(self, adata: AnnOrMuData) -> None:
        """
        Run a specified hyperparameter sweep for the associated model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been
            setup with the associated model class.
        """
