from __future__ import annotations

from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager


def get_minified_adata_scrna(
    adata_manager: AnnDataManager,
    keep_count_data: bool = False,
) -> AnnData:
    """Get a minified version of an :class:`~anndata.AnnData` or :class:`~mudata.MuData` object."""
    counts = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    return AnnData(
        X=counts if keep_count_data else csr_matrix(counts.shape),
        obs=adata_manager.adata.obs.copy(),
        var=adata_manager.adata.var.copy(),
        uns=adata_manager.adata.uns.copy(),
        obsm=adata_manager.adata.obsm.copy(),
        varm=adata_manager.adata.varm.copy(),
        obsp=adata_manager.adata.obsp.copy(),
        varp=adata_manager.adata.varp.copy(),
    )
