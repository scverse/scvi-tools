from anndata import AnnData
from mudata import MuData
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
) -> AnnData:
    """Returns a minified adata that works for most scrna models (such as SCVI, SCANVI).

    Parameters
    ----------
    adata
        Original adata, of which we to create a minified version.
    minified_data_type
        How to minify the data.
    """
    if minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
        raise NotImplementedError(f"Unknown MinifiedDataType: {minified_data_type}")

    all_zeros = csr_matrix(adata.X.shape)
    layers = {layer: all_zeros.copy() for layer in adata.layers}
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
