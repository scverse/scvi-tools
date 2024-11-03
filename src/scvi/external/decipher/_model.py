import logging
from collections.abc import Sequence

import numpy as np
import pyro
import torch
from anndata import AnnData

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField
from scvi.model.base import BaseModelClass, PyroSviTrainMixin
from scvi.train import PyroTrainingPlan
from scvi.utils import setup_anndata_dsp

from ._module import DecipherPyroModule
from ._trainingplan import DecipherTrainingPlan

logger = logging.getLogger(__name__)


class Decipher(PyroSviTrainMixin, BaseModelClass):
    """Decipher model for single-cell data analysis :cite:p:`Nazaret23`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via
        :meth:`~scvi.model.Decipher.setup_anndata`.
    dim_v
        Dimension of the interpretable latent space v.
    dim_z
        Dimension of the intermediate latent space z.
    layers_v_to_z
        Hidden layer sizes for the v to z decoder network.
    layers_z_to_x
        Hidden layer sizes for the z to x decoder network.
    beta
        Regularization parameter for the KL divergence.
    """

    _module_cls = DecipherPyroModule
    _training_plan_cls = DecipherTrainingPlan

    def __init__(self, adata: AnnData, **kwargs):
        pyro.clear_param_store()

        super().__init__(adata)

        dim_genes = self.summary_stats.n_vars

        self.module = self._module_cls(
            dim_genes,
            **kwargs,
        )

        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        **kwargs,
    ) -> AnnData | None:
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        early_stopping: bool = False,
        training_plan: PyroTrainingPlan | None = None,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        if "early_stopping_monitor" not in trainer_kwargs:
            trainer_kwargs["early_stopping_monitor"] = "nll_validation"
        datasplitter_kwargs = datasplitter_kwargs or {}
        if "drop_last" not in datasplitter_kwargs:
            datasplitter_kwargs["drop_last"] = True
        super().train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            device=device,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            training_plan=training_plan,
            datasplitter_kwargs=datasplitter_kwargs,
            **trainer_kwargs,
        )

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        give_z: bool = False,
    ) -> np.ndarray:
        """Get the latent representation of the data.

        Parameters
        ----------
        adata
            AnnData object with the data to get the latent representation of.
        indices
            Indices of the data to get the latent representation of.
        batch_size
            Batch size to use for the data loader.
        give_z
            Whether to return the intermediate latent space z or the top-level
            latent space v.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent_locs = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            x = torch.log1p(x)
            x = x.to(self.module.device)
            z_loc, _ = self.module.encoder_x_to_z(x)
            if give_z:
                latent_locs.append(z_loc)
            else:
                v_loc, _ = self.module.encoder_zx_to_v(torch.cat([z_loc, x], dim=-1))
                latent_locs.append(v_loc)
        return torch.cat(latent_locs).detach().cpu().numpy()
