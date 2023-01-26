import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from pytorch_lightning.callbacks import Callback

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
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.train import LoudEarlyStopping, TrainingPlan, TrainRunner
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)

B = 10


class CellAssign(UnsupervisedTrainingMixin, BaseModelClass):
    """
    Reimplementation of CellAssign for reference-based annotation :cite:p:`Zhang19`.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via :meth:`~scvi.external.CellAssign.setup_anndata`.
        The object should be subset to contain the same genes as the cell type marker dataframe.
    cell_type_markers
        Binary marker gene DataFrame of genes by cell types. Gene names corresponding to `adata.var_names`
        should be in DataFrame index, and cell type labels should be the columns.
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
    Size factors in the R implementation of CellAssign are computed using scran. An approximate approach
    computes the sum of UMI counts (library size) over all genes and divides by the mean library size.
    """

    def __init__(
        self,
        adata: AnnData,
        cell_type_markers: pd.DataFrame,
        **model_kwargs,
    ):
        try:
            cell_type_markers = cell_type_markers.loc[adata.var_names]
        except KeyError:
            raise KeyError(
                "Anndata and cell type markers do not contain the same genes."
            )
        super().__init__(adata)

        self.n_genes = self.summary_stats.n_vars
        self.cell_type_markers = cell_type_markers
        rho = torch.Tensor(cell_type_markers.to_numpy())
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(
                REGISTRY_KEYS.CAT_COVS_KEY
            ).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        adata = self._validate_anndata(adata)
        x = self.get_from_registry(adata, REGISTRY_KEYS.X_KEY)
        col_means = np.asarray(np.mean(x, 0)).ravel()  # (g)
        col_means_mu, col_means_std = np.mean(col_means), np.std(col_means)
        col_means_normalized = torch.Tensor((col_means - col_means_mu) / col_means_std)

        # compute basis means for phi - shape (B)
        basis_means = np.linspace(np.min(x), np.max(x), B)  # (B)

        self.module = CellAssignModule(
            n_genes=self.n_genes,
            rho=rho,
            basis_means=basis_means,
            b_g_0=col_means_normalized,
            n_batch=self.summary_stats.n_batch,
            n_cats_per_cov=n_cats_per_cov,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            **model_kwargs,
        )
        self._model_summary_string = (
            "CellAssign Model with params: \nn_genes: {}, n_labels: {}"
        ).format(
            self.n_genes,
            rho.shape[1],
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def predict(self) -> pd.DataFrame:
        """Predict soft cell type assignment probability for each cell."""
        adata = self._validate_anndata(None)
        scdl = self._make_data_loader(adata=adata)
        predictions = []
        for tensors in scdl:
            generative_inputs = self.module._get_generative_input(tensors, None)
            outputs = self.module.generative(**generative_inputs)
            gamma = outputs["gamma"]
            predictions += [gamma.cpu()]
        return pd.DataFrame(
            torch.cat(predictions).numpy(), columns=self.cell_type_markers.columns
        )

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 3e-3,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 1024,
        plan_kwargs: Optional[dict] = None,
        early_stopping: bool = True,
        early_stopping_patience: int = 15,
        early_stopping_min_delta: float = 0.0,
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
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.train.ClassifierTrainingPlan`. Keyword arguments passed to
        early_stopping
            Adds callback for early stopping on validation_loss
        early_stopping_patience
            Number of times early stopping metric can not improve over early_stopping_min_delta
        early_stopping_min_delta
            Threshold for counting an epoch torwards patience
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {"lr": lr, "weight_decay": 1e-10}
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        if "callbacks" in kwargs:
            kwargs["callbacks"] += [ClampCallback()]
        else:
            kwargs["callbacks"] = [ClampCallback()]

        if early_stopping:
            early_stopping_callback = [
                LoudEarlyStopping(
                    monitor="elbo_validation",
                    min_delta=early_stopping_min_delta,
                    patience=early_stopping_patience,
                    mode="min",
                )
            ]
            if "callbacks" in kwargs:
                kwargs["callbacks"] += early_stopping_callback
            else:
                kwargs["callbacks"] = early_stopping_callback
            kwargs["check_val_every_n_epoch"] = 1

        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = int(np.min([round((20000 / n_cells) * 400), 400]))

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = TrainingPlan(self.module, **plan_kwargs)
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **kwargs,
        )
        return runner()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        size_factor_key: str,
        batch_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
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
            CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
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
