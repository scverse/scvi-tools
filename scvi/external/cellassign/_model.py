import numpy as np
import pandas as pd
import torch

# import pdb
from anndata import AnnData

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
        sc_adata: AnnData,
        cell_type_markers: pd.DataFrame,
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super().__init__(sc_adata, use_gpu=use_gpu)
        # check that genes are the same in sc_adata and cell_type_markers
        if not sc_adata.var.index.equals(cell_type_markers.index):
            raise ValueError(
                "Genes must be the same in sc_adata and cell_type_markers (rho matrix)."
            )

        # reorder cell_type_markers according to order of genes in sc_adata.var.index
        cell_type_markers.reindex(sc_adata.var.index)

        self.n_genes = self.summary_stats["n_vars"]
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

    def predict(self, adata: AnnData) -> np.ndarray:
        """Predict soft cell type assignment probability for each cell."""
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata)
        latent = []
        for tensors in scdl:
            inference_inputs = self.model._get_inference_input(tensors)
            x, y = self.model.inference(**inference_inputs)
            outputs = self.model.generative(self, x, y)
            gamma = outputs["gamma"]
            latent += [gamma.cpu()]
        return np.array(torch.cat(latent))

    @property
    def _task_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    @property
    def _plan_class(self):
        return TrainingPlan
