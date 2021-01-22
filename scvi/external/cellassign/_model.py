import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi.dataloaders import ScviDataLoader
from scvi.external.cellassign._module import CellAssignModule
from scvi.lightning import VAETask
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
        sc_adata: AnnData,
        cell_type_markers: pd.DataFrame,
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super().__init__(sc_adata, use_gpu=use_gpu)
        self.n_genes = self.summary_stats["n_vars"]
        self.n_labels = self.summary_stats["n_labels"]

        # might need to reorganize the df to correspond to the
        # correct order of labels
        rho = torch.Tensor(cell_type_markers.to_numpy())
        self.model = CellAssignModule(
            n_genes=self.n_genes,
            n_labels=self.n_labels,
            rho=rho,
            **model_kwargs,
        )
        self._model_summary_string = (
            "CellAssign Model with params: \nn_genes: {}, n_labels: {}"
        ).format(
            self.n_genes,
            self.n_labels,
        )
        self.init_params_ = self._get_init_params(locals())

    def predict(self) -> np.ndarray:
        """Predict soft cell type assignment probability for each cell."""
        raise NotImplementedError

    @property
    def _task_class(self):
        return VAETask

    @property
    def _data_loader_cls(self):
        return ScviDataLoader
