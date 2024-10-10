from collections.abc import Sequence
from typing import TYPE_CHECKING

import logging

from anndata import AnnData
import pyro

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField, CategoricalJointObsField
from scvi.utils import setup_anndata_dsp
from scvi.train import PyroTrainingPlan

from scvi.model.base import BaseModelClass, PyroSviTrainMixin

from ._module import DecipherPyroModule

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

logger = logging.getLogger(__name__)


class Decipher(PyroSviTrainMixin, BaseModelClass):
    _module_cls = DecipherPyroModule

    def __init__(self, adata: AnnData, **kwargs):
        pyro.clear_param_store()

        super().__init__(adata)

        dim_genes = self.summary_stats.n_vars

        self.module = self._module_cls(
            dim_genes,
            **kwargs,
        )

        self.init_params = self._get_init_params(locals())

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
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
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
        lr: float | None = None,
        training_plan: PyroTrainingPlan | None = None,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        optim_kwargs = trainer_kwargs.pop("optim_kwargs", {})
        optim_kwargs.update({"lr": lr or 5e-3, "weight_decay": 1e-4})
        optim = pyro.optim.ClippedAdam(optim_kwargs)
        plan_kwargs = plan_kwargs or {}
        plan_kwargs.update({"optim": optim})
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
