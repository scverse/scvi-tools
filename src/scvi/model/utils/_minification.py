from __future__ import annotations

from typing import TYPE_CHECKING

from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS
from scvi.data._constants import (
    _ADATA_MINIFY_TYPE_UNS_KEY,
    _SCVI_UUID_KEY,
    ADATA_MINIFY_TYPE,
)

if TYPE_CHECKING:
    from mudata import MuData

    from scvi._types import MinifiedDataType
    from scvi.data import AnnDataManager


def get_minified_adata_scrna(
    adata_manager: AnnDataManager,
    keep_count_data: bool = False,
) -> AnnData:
    """Get a minified version of an :class:`~anndata.AnnData` or :class:`~mudata.MuData` object."""
    counts = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    all_zeros = csr_matrix(counts.shape)
    return AnnData(
        X=counts if keep_count_data else all_zeros,
        layers={layer: all_zeros.copy() for layer in adata_manager.adata.layers},
        obs=adata_manager.adata.obs.copy(),
        var=adata_manager.adata.var.copy(),
        uns=adata_manager.adata.uns.copy(),
        obsm=adata_manager.adata.obsm.copy(),
        varm=adata_manager.adata.varm.copy(),
        obsp=adata_manager.adata.obsp.copy(),
        varp=adata_manager.adata.varp.copy(),
    )


def get_minified_mudata(
    mdata: MuData,
    minified_data_type: MinifiedDataType,
) -> MuData:
    """Returns a minified adata that works for most multi modality models (MULTIVI, TOTALVI).

    Parameters
    ----------
    mdata
        Original adata, of which we to create a minified version.
    minified_data_type
        How to minify the data.
    """
    if minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
        raise NotImplementedError(f"Unknown MinifiedDataType: {minified_data_type}")

    bdata = mdata.copy()
    for modality in mdata.mod_names:
        all_zeros = csr_matrix(mdata[modality].X.shape)
        bdata[modality].X = all_zeros
        if len(mdata[modality].layers) > 0:
            layers = {layer: all_zeros for layer in mdata[modality].layers}
            bdata[modality].layers = layers
    # Remove scvi uuid key to make bdata fresh w.r.t. the model's manager
    del bdata.uns[_SCVI_UUID_KEY]
    bdata.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = minified_data_type
    return bdata
