import logging

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

import scvi
from scvi import _CONSTANTS
from scvi.data import register_tensor_from_anndata
from scvi.external.cellassign._module import CellAssignModule
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin

logger = logging.getLogger(__name__)

B = 10


class CellAssign(UnsupervisedTrainingMixin, BaseModelClass):
    """
    Reimplementation of CellAssign for reference-based annotation [Zhang19]_.

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
    >>> adata = scvi.data.read_h5ad(path_to_anndata)
    >>> marker_gene_mat = pd.read_csv(path_to_marker_gene_csv)
    >>> bdata = adata[:, adata.var.index.isin(marker_gene_mat.index)].copy()
    >>> scvi.data.setup_anndata(bdata)
    >>> model = CellAssign(bdata, marker_gene_mat,size_factor_key='S')
    >>> model.train(max_epochs=100, batch_size=1024)
    >>> predictions = model.predict(bdata)
    """

    def __init__(
        self,
        adata: AnnData,
        cell_type_markers: pd.DataFrame,
        size_factor_key: str,
        **model_kwargs,
    ):
        super().__init__(adata)
        # check that genes are the same in cell_type_markers are present in the anndata
        # anndata may have more
        adata.var = adata.var.sort_index()
        cell_type_markers = cell_type_markers.sort_index()
        if not adata.var.index.equals(cell_type_markers.index):
            raise ValueError(
                "Anndata and cell type markers do not contain the same genes."
            )
        cell_type_markers = cell_type_markers.loc[adata.var_names]

        register_tensor_from_anndata(adata, "_size_factor", "obs", size_factor_key)

        self.n_genes = self.summary_stats["n_vars"]
        self.cell_type_markers = cell_type_markers
        rho = torch.Tensor(cell_type_markers.to_numpy())
        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )

        x = scvi.data.get_from_registry(adata, _CONSTANTS.X_KEY)
        col_means = np.mean(x, 0)  # (g)
        col_means_mu, col_means_std = np.mean(col_means), np.std(col_means)
        col_means_normalized = torch.Tensor((col_means - col_means_mu) / col_means_std)

        # compute basis means for phi - shape (B)
        basis_means = np.linspace(np.min(x), np.max(x), B)  # (B)

        self.module = CellAssignModule(
            n_genes=self.n_genes,
            rho=rho,
            b_g_0=col_means_normalized,
            basis_means=basis_means,
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
    def predict(self) -> np.ndarray:
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
            np.array(torch.cat(predictions)), columns=self.cell_type_markers.columns
        )
