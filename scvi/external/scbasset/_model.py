import logging
from typing import Literal, Optional, Union

import numpy as np
import torch
from anndata import AnnData

from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalVarField, LayerField, ObsmField
from scvi.dataloaders import DataSplitter
from scvi.external.scbasset._module import REGISTRY_KEYS, ScBassetModule
from scvi.model.base import BaseModelClass
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)


class SCBASSET(BaseModelClass):
    """
    Reimplementation of ScBasset :cite:p:`Yuan2022` for representation learning of scATAC-seq data.

    This implementation is EXPERIMENTAL. We are working to measure the performance of this model
    compared to the original.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via :meth:`~scvi.external.SCBASSET.setup_anndata`.
    **model_kwargs
        Keyword args for :class:`~scvi.external.scbasset.ScBassetModule`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_sc_anndata)
    >>> scvi.data.add_dna_sequence(adata)
    >>> adata = adata.transpose() # regions by cells
    >>> scvi.external.SCBASSET.setup_anndata(adata, dna_code_key="dna_code")
    >>> model = scvi.external.SCBASSET(adata)
    >>> model.train()
    >>> adata.varm["X_scbasset"] = model.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        **model_kwargs,
    ):
        super().__init__(adata)
        self.n_cells = self.summary_stats.n_vars
        self.n_regions = adata.n_obs
        self.n_batch = self.summary_stats.n_batch
        batch_ids = self.adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY)
        self.module = ScBassetModule(
            n_cells=self.n_cells,
            batch_ids=torch.tensor(batch_ids).long() if batch_ids.sum() > 0 else None,
            **model_kwargs,
        )
        self._model_summary_string = (
            "ScBasset Model with params: \nn_regions: {}, n_batch: {}, n_cells: {}"
        ).format(
            self.n_regions,
            self.n_batch,
            self.n_cells,
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 1000,
        lr: float = 0.01,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        early_stopping: bool = True,
        early_stopping_monitor: str = "auroc_train",
        early_stopping_mode: Literal["min", "max"] = "max",
        early_stopping_min_delta: float = 1e-6,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        early_stopping_monitor
            Metric logged during validation set epoch. The available metrics will depend on
            the training plan class used. We list the most common options here in the typing.
        early_stopping_mode
            In 'min' mode, training will stop when the quantity monitored has stopped decreasing
            and in 'max' mode it will stop when the quantity monitored has stopped increasing.
        early_stopping_min_delta
            Minimum change in the monitored quantity to qualify as an improvement,
            i.e. an absolute change of less than min_delta, will count as no improvement.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        custom_plan_kwargs = dict(
            optimizer="Custom",
            optimizer_creator=lambda p: torch.optim.Adam(
                p, lr=lr, betas=(0.95, 0.9995)
            ),
        )
        if plan_kwargs is not None:
            custom_plan_kwargs.update(plan_kwargs)

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
            # We don't want to dataload the batch ids into the module
            data_and_attributes={
                REGISTRY_KEYS.X_KEY: np.float32,
                REGISTRY_KEYS.DNA_CODE_KEY: np.int64,
            },
        )
        training_plan = TrainingPlan(self.module, **custom_plan_kwargs)

        es = {
            "early_stopping": early_stopping,
            "early_stopping_monitor": early_stopping_monitor,
            "early_stopping_mode": early_stopping_mode,
            "early_stopping_min_delta": early_stopping_min_delta,
        }
        for k, v in es.items():
            trainer_kwargs[k] = (
                v if k not in trainer_kwargs.keys() else trainer_kwargs[k]
            )
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()

    @torch.inference_mode()
    def get_latent_representation(self) -> np.ndarray:
        """
        Returns the latent representation of the cells.

        Returns
        -------
        latent representation (n_cells, n_latent)
        """
        return self.module.cell_embedding.cpu().numpy().T

    @torch.inference_mode()
    def get_cell_bias(self) -> np.ndarray:
        """
        Returns the cell-specific bias term.

        Returns
        -------
        bias (n_cells,)
        """
        return self.module.cell_bias.cpu().numpy()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        dna_code_key: str,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        dna_code_key
            Key in `adata.obsm` with dna sequences encoded as integer code.
        %(param_layer)s
        batch_key
            key in `adata.var` for batch information. Categories will automatically be converted into integer
            categories and saved to `adata.var['_scvi_batch']`. If `None`, assigns the same batch to all the data.

        Notes
        -----
        The adata object should be in the regions by cells format. This is due to scBasset
        considering regions as observations and cells as variables. This can be simply achieved
        by transposing the data, `bdata = adata.transpose()`.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            ObsmField(REGISTRY_KEYS.DNA_CODE_KEY, dna_code_key, is_count_data=True),
            CategoricalVarField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
