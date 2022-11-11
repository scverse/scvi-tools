from scvi.model.base import BaseModelClass


class TunerManager:
    """
    Internal manager for validation of inputs from :class:`~scvi.autotune.ModelTuner`.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune
        hyperparameters. See :class:`~scvi.autotune.ModelTuner` for
        supported model classes.
    """

    def __init__(self, model_cls: BaseModelClass):
        pass
