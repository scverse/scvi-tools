import io
import logging
import warnings
from contextlib import redirect_stdout
from typing import Optional

import numpy as np
import torch
from anndata import AnnData

from scvi import _CONSTANTS
from scvi.compose import auto_move_data
from scvi.data import get_from_registry, setup_anndata
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import ClassifierTrainingPlan
from scvi.model import SCVI
from scvi.model.base import BaseModelClass
from scvi.modules import Classifier

logger = logging.getLogger(__name__)

LABELS_KEY = "_solo_doub_sim"


class SOLO(BaseModelClass):
    """
    Doublet detection in scRNA-seq [Bernstein19]_.

    Most users will initialize the model using the class method
    :func:`~scvi.external.SOLO.from_scvi_model`, which takes as
    input a pre-trained :class:`~scvi.model.SCVI` object.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
        Object should contain latent representation of real cells and doublets as `adata.X`.
        Object should also be registered, using `.X` and `labels_key="_solo_doub_sim"`.
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
        use_gpu: bool = True,
        **model_kwargs,
    ):
        # TODO, catch user warning here and logger warning
        # about non count data
        super().__init__(adata, use_gpu=use_gpu)

        self.module = Classifier(
            n_input=self.summary_stats["n_vars"],
            n_labels=2,
            logits=True,
            **model_kwargs,
        )
        self._model_summary_string = "Solo model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_scvi_model(cls, scvi_model: SCVI, adata: Optional[AnnData] = None):
        """
        Instantiate a SOLO model from an scvi model.

        Parameters
        ----------
        scvi_model
            Pre-trained model of :class:`~scvi.model.SCVI`. This model
            should have been trained on data comprising one lane. The
            adata object used to initialize this model should have only
            been setup with count data, i.e., no `batch_key`,
            `labels_key`, etc.
        adata
            Optional anndata to use that is compatible with scvi_model.

        Returns
        -------
        SOLO model
        """
        _validate_scvi_model(scvi_model)
        doublet_adata = cls.create_doublets(scvi_model.adata)

        # if model is using observed lib size, needs to get lib sample
        # which is just observed lib size on log scale
        give_mean_lib = not scvi_model.module.use_observed_lib_size

        # get latent representations and make input anndata
        latent_rep = scvi_model.get_latent_representation()
        lib_size = scvi_model.get_latent_library_size(give_mean=give_mean_lib)
        latent_adata = AnnData(np.concatenate([latent_rep, lib_size], axis=1))
        latent_adata.obs[LABELS_KEY] = "singlet"

        logger.info("Creating doublets, preparing SOLO model.")
        f = io.StringIO()
        with redirect_stdout(f):
            setup_anndata(doublet_adata)
            doublet_latent_rep = scvi_model.get_latent_representation(doublet_adata)
            doublet_lib_size = scvi_model.get_latent_library_size(
                doublet_adata, give_mean=give_mean_lib
            )
            doublet_adata = AnnData(
                np.concatenate([doublet_latent_rep, doublet_lib_size], axis=1)
            )
            doublet_adata.obs[LABELS_KEY] = "doublet"

            full_adata = latent_adata.concatenate(doublet_adata)
            setup_anndata(full_adata, labels_key=LABELS_KEY)
        return cls(full_adata)

    @staticmethod
    def create_doublets(
        adata: AnnData, seed: int = 1, doublet_ratio: int = 2
    ) -> AnnData:
        """Simulate doublets."""
        num_doublets = doublet_ratio * adata.n_obs

        # counts can be in many locations, this uses where it was registered in setup
        x = get_from_registry(adata, _CONSTANTS.X_KEY)

        # TODO: needs a random state so it's reproducible
        parent_inds = np.random.choice(adata.n_obs, size=(num_doublets, 2))
        doublets = x[parent_inds[:, 0]] + x[parent_inds[:, 1]]

        doublets_ad = AnnData(doublets)
        doublets_ad.var_names = adata.var_names
        doublets_ad.obs_names = [
            "sim_doublet_{}".format(i) for i in range(num_doublets)
        ]

        return doublets_ad

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

    def _set_training_plan(self, plan_class, plan_kwargs):
        """Set the _training_plan attribute."""
        if plan_class is None:
            plan_class = self._plan_class

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        self._training_plan = plan_class(self.module, **plan_kwargs)

    @torch.no_grad()
    def predict(self, soft: bool = True):
        """
        Return doublet predictions.

        Parameters
        ----------
        soft
            Return probabilities instead of class label
        """
        adata = self._validate_anndata(None)

        scdl = self._make_data_loader(
            adata=adata,
        )

        @auto_move_data
        def auto_forward(module, x):
            return module(x)

        y_pred = []
        for _, tensors in enumerate(scdl):
            x = tensors[_CONSTANTS.X_KEY]
            pred = auto_forward(self.module, x)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.cpu())

        y_pred = torch.cat(y_pred).numpy()

        # TODO: make it a dataframe with nice columns

        return y_pred

    @property
    def _plan_class(self):
        return ClassifierTrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader


def _validate_scvi_model(scvi_model: SCVI):
    if scvi_model.summary_stats["n_batch"] > 1:
        warnings.warn(
            "The SCVI model should only be trained on one lane of data. Performance may suffer.",
            UserWarning,
        )
