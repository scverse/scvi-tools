import io
from contextlib import redirect_stdout
from typing import Optional

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse

from scvi import _CONSTANTS
from scvi.data import get_from_registry, setup_anndata
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import ClassifierTrainingPlan
from scvi.model import SCVI
from scvi.model.base import BaseModelClass
from scvi.modules import Classifier


LABELS_KEY = "_solo_doub_sim"


class SOLO(BaseModelClass):
    """
    Reimplementation of Solo for doublet detection.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    latent_obsm_key
        Key in adata.obsm containing SCVI latent representation
    use_gpu
        Use the GPU or not.
    **classifier_kwargs
        Keyword args for :class:`~scvi.modules.Classifier`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.data.setup_anndata(adata)
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> solo = scvi.external.SOLO.from_scvi_model(vae)
    >>> solo.train()
    >>> solo.predict()
    """

    def __init__(
        self,
        adata: AnnData,
        latent_obsm_key: str,
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super().__init__(adata, use_gpu=use_gpu)

        latent_rep = adata.obsm[latent_obsm_key]
        x = get_from_registry(adata, _CONSTANTS.X_KEY)
        if issparse(x):
            library = np.log(x.sum(1).A.ravel().reshape(-1, 1))
        else:
            library = np.log(x.sum(1).reshape(-1, 1))

        self._classifier_data = AnnData(np.concatenate([latent_rep, x], axis=1))
        self._classifier_data = self.add_doublets(self._classifier_data)

        f = io.StringIO()
        with redirect_stdout(f):
            setup_anndata(self._classifier_data, labels_key=LABELS_KEY)

        self.module = Classifier(
            n_input=self._scvi_model.n_latent + 1,
            n_labels=2,
            logits=True,
            **model_kwargs,
        )
        self._model_summary_string = "Solo model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_scvi_model(cls, scvi_model: SCVI):
        """Instantiate a model from an scvi model."""
        latent_obsm_key = "_solo_latent_rep"
        _validate_scvi_model(scvi_model)
        latent_rep = scvi_model.get_latent_representation()
        scvi_model.adata.obsm[latent_obsm_key] = latent_rep
        return cls(scvi_model.adata, latent_obsm_key)

    @staticmethod
    def add_doublets(adata: AnnData):
        """Simulate and add doublets to anndata."""

        # TODO needs a random seed, use scvi.settings.seed
        # add an obs column with name LABELS_KEY of strings, "singlet", "doublet"
        return adata

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 1e-3,
        use_gpu: Optional[bool] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.ClassifierTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        update_dict = {
            "lr": lr,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    def predict(self):
        raise NotImplementedError

    @property
    def _plan_class(self):
        return ClassifierTrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader


def _validate_scvi_model(scvi_model: SCVI):
    # check that it was trained on one batch
    pass
