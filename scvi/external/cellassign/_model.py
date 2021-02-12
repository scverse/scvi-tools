import numpy as np
import pandas as pd
import torch

# import pdb
from anndata import AnnData
from scipy.sparse import issparse

from scvi import _CONSTANTS
from scvi.data import get_from_registry, register_tensor_from_anndata
from scvi.dataloaders import AnnDataLoader
from scvi.external.cellassign._module import CellAssignModule
from scvi.lightning import TrainingPlan
from scvi.model.base import BaseModelClass


class CellAssign(BaseModelClass):
    """
    Reimplementation of CellAssign for reference-based annotation.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    cell_type_markers
        Binary marker gene matrix
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.CellAssignModels`

    Examples
    --------
    >>> # TODO
    """

    def __init__(
        self,
        adata: AnnData,
        cell_type_markers: pd.DataFrame,
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super().__init__(adata, use_gpu=use_gpu)
        # check that genes are the same in adata and cell_type_markers
        if not adata.var.index.sort_values().equals(
            cell_type_markers.index.sort_values()
        ):
            raise ValueError(
                "Genes must be the same in adata and cell_type_markers (rho matrix)."
            )

        counts = get_from_registry(adata, _CONSTANTS.X_KEY)
        if issparse(counts):
            size_factor = counts.sum(1).A.ravel()
        else:
            size_factor = counts.sum(1).ravel()
        adata.obs["_size_factor"] = size_factor
        register_tensor_from_anndata(adata, "_size_factor", "obs", "_size_factor")

        self.n_genes = self.summary_stats["n_vars"]
        self.cell_type_markers = cell_type_markers
        rho = torch.Tensor(cell_type_markers.to_numpy())
        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )

        self.module = CellAssignModule(
            n_genes=self.n_genes,
            rho=rho,
            n_batch=self.summary_stats["n_batch"],
            n_cats_per_cov=n_cats_per_cov,
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            **model_kwargs,
        )
        self._model_summary_string = (
            "CellAssign Model with params: \nn_genes: {}, n_labels: {}"
        ).format(
            self.n_genes,
            rho.shape[1],
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.no_grad()
    def predict(self, adata: AnnData) -> np.ndarray:
        """Predict soft cell type assignment probability for each cell."""
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata)
        predictions = []
        for tensors in scdl:
            generative_inputs = self.module._get_generative_input(tensors, None)
            outputs = self.module.generative(**generative_inputs)
            gamma = outputs["gamma"]
            predictions += [gamma.cpu()]
        return pd.DataFrame(
            np.array(torch.cat(predictions)), columns=self.cell_type_markers.columns
        )

    @property
    def _task_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    @property
    def _plan_class(self):
        return TrainingPlan
