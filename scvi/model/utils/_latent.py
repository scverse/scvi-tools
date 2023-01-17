from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi._types import LatentDataType
from scvi.data._constants import _ADATA_LATENT_UNS_KEY, _SCVI_UUID_KEY


def get_reduced_adata(
    adata: AnnData,
    mode: LatentDataType,
) -> AnnData:
    """Return a minimal anndata object with the latent representation."""
    all_zeros = csr_matrix(adata.X.shape)
    layers = {layer: all_zeros.copy() for layer in adata.layers}
    bdata = AnnData(
        X=all_zeros,
        layers=layers,
        uns=adata.uns,
        obs=adata.obs,
        var=adata.var,
        varm=adata.varm,
        obsm=adata.obsm,
        obsp=adata.obsp,
    )
    # Remove scvi uuid key to make bdata fresh w.r.t. the model's manager
    del bdata.uns[_SCVI_UUID_KEY]
    bdata.uns[_ADATA_LATENT_UNS_KEY] = mode
    return bdata
