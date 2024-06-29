from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi._types import MinifiedDataType
from scvi.data._constants import (
    _ADATA_MINIFY_TYPE_UNS_KEY,
    _SCVI_UUID_KEY,
    ADATA_MINIFY_TYPE,
)


def get_minified_adata_scrna(
    adata: AnnData,
    minified_data_type: MinifiedDataType,
    keep_count_data: bool = False,
) -> AnnData:
    """Returns a minified adata that works for most scrna models (such as SCVI, SCANVI).

    Parameters
    ----------
    adata
        Original adata, of which we to create a minified version.
    minified_data_type
        How to minify the data.
    keep_count_data
        Whether to keep the count data and only store additionally the latent posterior.
    """
    if not ADATA_MINIFY_TYPE.__contains__(minified_data_type):
        raise NotImplementedError(f"Unknown MinifiedDataType: {minified_data_type}")

    if not keep_count_data:
        all_zeros = csr_matrix(adata.X.shape)
        layers = {layer: all_zeros.copy() for layer in adata.layers}
    else:
        all_zeros = adata.X
        layers = adata.layers
    bdata = AnnData(
        X=all_zeros,
        layers=layers,
        uns=adata.uns.copy(),
        obs=adata.obs,
        var=adata.var,
        varm=adata.varm,
        obsm=adata.obsm,
        obsp=adata.obsp,
    )
    # Remove scvi uuid key to make bdata fresh w.r.t. the model's manager
    del bdata.uns[_SCVI_UUID_KEY]
    bdata.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = minified_data_type
    return bdata
