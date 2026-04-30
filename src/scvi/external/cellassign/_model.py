from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from lightning.pytorch.callbacks import Callback

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.dataloaders import DataSplitter
from scvi.external.cellassign._module import CellAssignModule
from scvi.model._utils import get_max_epochs_heuristic
from scvi.model.base import BaseModelClass, RNASeqMixin, UnsupervisedTrainingMixin
from scvi.model.base._training_mixin import _set_datamodule_split
from scvi.train import LoudEarlyStopping, TrainingPlan, TrainRunner
from scvi.train._config import merge_kwargs
from scvi.utils import dependencies, setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

B = 10


class CellAssign(UnsupervisedTrainingMixin, RNASeqMixin, BaseModelClass):
    """Reimplementation of CellAssign for reference-based annotation :cite:p:`Zhang19`.

    Original implementation: https://github.com/irrationone/cellassign.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via
        :meth:`~scvi.external.CellAssign.setup_anndata`. The object should be subset to contain the
        same genes as the cell type marker dataframe.
    cell_type_markers
        Binary marker gene DataFrame of genes by cell types. Gene names corresponding to
        `adata.var_names` should be in the DataFrame index,
        and cell type labels should be the columns.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cellassign.CellAssignModule`

    Examples
    --------
    >>> adata = scvi.data.read_h5ad(path_to_anndata)
    >>> library_size = adata.X.sum(1)
    >>> adata.obs["size_factor"] = library_size / np.mean(library_size)
    >>> marker_gene_mat = pd.read_csv(path_to_marker_gene_csv)
    >>> bdata = adata[:, adata.var.index.isin(marker_gene_mat.index)].copy()
    >>> CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
    >>> model = CellAssign(bdata, marker_gene_mat)
    >>> model.train()
    >>> predictions = model.predict(bdata)

    Notes
    -----
    Size factors in the R implementation of CellAssign are computed using scran. An approximate
    approach computes the sum of UMI counts (library size) over all genes and divides by the mean
    library size.

    See further usage examples in the following tutorial:

    1. :doc:`/tutorials/notebooks/scrna/cellassign_tutorial`
    """

    def __init__(
        self,
        adata: AnnData | None = None,
        cell_type_markers: pd.DataFrame | None = None,
        registry: dict | None = None,
        col_means: np.ndarray | None = None,
        basis_means: np.ndarray | None = None,
        **model_kwargs,
    ):
        if adata is not None:
            try:
                cell_type_markers = cell_type_markers.loc[adata.var_names]
            except KeyError as err:
                raise KeyError(
                    "Anndata and cell type markers do not contain the same genes."
                ) from err

        assert not cell_type_markers.index.has_duplicates, (
            "There are duplicates in cell type markers (rows in cell_type_markers)"
        )

        super().__init__(adata, registry)

        self.n_genes = self.summary_stats.n_vars
        self.cell_type_markers = cell_type_markers
        rho = torch.Tensor(cell_type_markers.to_numpy())

        if adata is not None:
            n_cats_per_cov = (
                self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
                if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
                else None
            )
            adata = self._validate_anndata(adata)
            x = self.get_from_registry(adata, REGISTRY_KEYS.X_KEY)
            col_means_arr = np.asarray(np.mean(x, 0)).ravel()
            basis_means_arr = np.linspace(np.min(x), np.max(x), B)
        else:
            if col_means is None or basis_means is None:
                raise ValueError(
                    "col_means and basis_means must be provided when using the registry path. "
                    "Use CellAssign.setup_annbatch(...) which computes them automatically."
                )
            cat_cov_sr = self.registry["field_registries"]["extra_categorical_covs"][
                "state_registry"
            ]
            n_cats_per_cov = tuple(cat_cov_sr["n_cats_per_key"]) if cat_cov_sr else None
            col_means_arr = np.asarray(col_means).ravel()
            basis_means_arr = np.asarray(basis_means)

        col_means_mu, col_means_std = np.mean(col_means_arr), np.std(col_means_arr)
        col_means_normalized = torch.Tensor((col_means_arr - col_means_mu) / col_means_std)

        self.module = CellAssignModule(
            n_genes=self.n_genes,
            rho=rho,
            basis_means=basis_means_arr,
            b_g_0=col_means_normalized,
            n_batch=self.summary_stats.n_batch,
            n_cats_per_cov=n_cats_per_cov,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            **model_kwargs,
        )
        self._model_summary_string = (
            f"CellAssign Model with params: \nn_genes: {self.n_genes}, n_labels: {rho.shape[1]}"
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def predict(self, dataloader=None) -> pd.DataFrame:
        """Predict soft cell type assignment probability for each cell."""
        if dataloader is None:
            adata = self._validate_anndata(None)
            scdl = self._make_data_loader(adata=adata)
        else:
            scdl = dataloader
        predictions = []
        for tensors in scdl:
            generative_inputs = self.module._get_generative_input(tensors, None)
            outputs = self.module.generative(**generative_inputs)
            gamma = outputs["gamma"]
            predictions += [gamma.cpu()]
        return pd.DataFrame(torch.cat(predictions).numpy(), columns=self.cell_type_markers.columns)

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 400,
        lr: float = 3e-3,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 1024,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        early_stopping: bool = True,
        early_stopping_patience: int = 15,
        early_stopping_warmup_epochs: int = 0,
        early_stopping_min_delta: float = 0.0,
        datamodule=None,
        **kwargs,
    ):
        """Trains the model.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of the training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`.
        early_stopping
            Adds callback for early stopping on validation_loss
        early_stopping_patience
            Number of times early stopping metric cannot improve over early_stopping_min_delta
        early_stopping_warmup_epochs
            Wait for a certain number of warm-up epochs before the early stopping starts monitoring
        early_stopping_min_delta
            Threshold for counting an epoch towards patience
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        datamodule
            Optional AnnbatchDataModule for streaming training (from setup_annbatch).
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {"lr": lr, "weight_decay": 1e-10}
        plan_kwargs = merge_kwargs(None, plan_kwargs, name="plan")
        plan_kwargs.update(update_dict)

        if "callbacks" in kwargs:
            kwargs["callbacks"] += [ClampCallback()]
        else:
            kwargs["callbacks"] = [ClampCallback()]

        if datamodule is not None:
            # Annbatch streaming path: use datamodule directly
            _set_datamodule_split(datamodule, train_size, validation_size, shuffle_set_split)
            if max_epochs is None:
                max_epochs = get_max_epochs_heuristic(datamodule.n_obs)
            training_plan = TrainingPlan(self.module, **plan_kwargs)
            runner = TrainRunner(
                self,
                training_plan=training_plan,
                data_splitter=datamodule,
                max_epochs=max_epochs,
                accelerator=accelerator,
                devices=devices,
                **kwargs,
            )
            return runner()

        # Standard AnnData path
        datasplitter_kwargs = datasplitter_kwargs or {}

        if early_stopping:
            early_stopping_callback = [
                LoudEarlyStopping(
                    monitor="elbo_validation",
                    min_delta=early_stopping_min_delta,
                    patience=early_stopping_patience,
                    mode="min",
                    warmup_epochs=early_stopping_warmup_epochs,
                )
            ]
            if "callbacks" in kwargs:
                kwargs["callbacks"] += early_stopping_callback
            else:
                kwargs["callbacks"] = early_stopping_callback
            kwargs["check_val_every_n_epoch"] = 1

        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            shuffle_set_split=shuffle_set_split,
            **datasplitter_kwargs,
        )
        training_plan = TrainingPlan(self.module, **plan_kwargs)
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            **kwargs,
        )
        return runner()

    @classmethod
    @dependencies("annbatch", "zarr")
    def setup_annbatch(
        cls,
        collection_path: str,
        paths: list[str] | None = None,
        cell_type_markers=None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        layer: str | None = None,
        batch_size: int = 4096,
        chunk_size: int = 256,
        preload_nchunks: int = 32,
        preload_to_gpu: bool = False,
        dataset_size: int = 2_097_152,
        shuffle: bool = False,
        var_subset: list[str] | None = None,
        rebuild: bool = False,
        **kwargs,
    ):
        """Build annbatch collection and compute global statistics for CellAssign init.

        Returns
        -------
        tuple of (AnnbatchDataModule, col_means: np.ndarray, basis_means: np.ndarray)
        """
        import numpy as np

        dm = super().setup_annbatch(
            collection_path=collection_path,
            paths=paths,
            batch_key=batch_key,
            labels_key=labels_key,
            layer=layer,
            batch_size=batch_size,
            chunk_size=chunk_size,
            preload_nchunks=preload_nchunks,
            preload_to_gpu=preload_to_gpu,
            dataset_size=dataset_size,
            shuffle=shuffle,
            var_subset=var_subset,
            rebuild=rebuild,
            **kwargs,
        )

        # Stream through all data to compute per-gene mean, global min, global max
        col_sums = np.zeros(dm.n_vars, dtype=np.float64)
        global_min = np.inf
        global_max = -np.inf
        n_cells = 0

        for batch in dm.inference_dataloader():
            x = batch["X"].numpy()  # (batch_size, n_vars)
            col_sums += x.sum(axis=0)
            global_min = min(global_min, float(x.min()))
            global_max = max(global_max, float(x.max()))
            n_cells += x.shape[0]

        col_means = (col_sums / n_cells).astype(np.float32)
        basis_means = np.linspace(global_min, global_max, B).astype(np.float32)

        return dm, col_means, basis_means

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        size_factor_key: str,
        batch_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        layer: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        size_factor_key
            key in `adata.obs` with continuous valued size factors.
        %(param_batch_key)s
        %(param_layer)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


class ClampCallback(Callback):
    """Clamp callback."""

    def __init__(self):
        super().__init__()

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        """Clamp parameters."""
        with torch.inference_mode():
            pl_module.module.delta_log.clamp_(np.log(pl_module.module.min_delta))
        super().on_train_batch_end(trainer, pl_module, outputs, batch, batch_idx)
