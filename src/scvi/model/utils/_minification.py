from __future__ import annotations

from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS
from scvi._types import MinifiedDataType
from scvi.data import AnnDataManager
from scvi.data._constants import (
    _ADATA_MINIFY_TYPE_UNS_KEY,
    _SCVI_UUID_KEY,
    ADATA_MINIFY_TYPE,
)


def get_minified_adata_scrna(
    adata_manager: AnnDataManager,
    minified_data_type: MinifiedDataType = ADATA_MINIFY_TYPE.LATENT_POSTERIOR,
    keep_count_data: bool = False,
) -> AnnData:
    """Get a minified version of an :class:`~anndata.AnnData` or :class:`~mudata.MuData` object."""
    if minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
        raise NotImplementedError(f"Minification method {minified_data_type} is not supported.")

    counts = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    mini_adata = AnnData(
        X=counts if keep_count_data else csr_matrix(counts.shape),
        obs=adata_manager.adata.obs.copy(),
        var=adata_manager.adata.var.copy(),
        uns=adata_manager.adata.uns.copy(),
        obsm=adata_manager.adata.obsm.copy(),
        varm=adata_manager.adata.varm.copy(),
        obsp=adata_manager.adata.obsp.copy(),
        varp=adata_manager.adata.varp.copy(),
    )
    del mini_adata.uns[_SCVI_UUID_KEY]
    mini_adata.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = minified_data_type

    return mini_adata
