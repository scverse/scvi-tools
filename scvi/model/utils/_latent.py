from anndata import AnnData
from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS
from scvi._types import LatentDataType
from scvi.data._constants import _ADATA_LATENT_UNS_KEY
from scvi.data.fields import ObsmField, StringUnsField


def scvi_get_latent_adata_from_adata(
    adata: AnnData,
    mode: LatentDataType,
    scvi_latent_qzm_key: str,
    scvi_latent_qzv_key: str,
    use_latent_qzm_key: str = "X_latent_qzm",
    use_latent_qzv_key: str = "X_latent_qzv",
):
    """TODO add docstring"""
    if mode == "dist":
        adata.obsm[scvi_latent_qzm_key] = adata.obsm[use_latent_qzm_key]
        adata.obsm[scvi_latent_qzv_key] = adata.obsm[use_latent_qzv_key]
    else:
        raise ValueError(f"Unknown latent mode: {mode}")
    adata.uns[_ADATA_LATENT_UNS_KEY] = mode
    del adata.raw
    all_zeros = csr_matrix(adata.X.shape)
    adata.X = all_zeros.copy()
    adata.layers = {layer: all_zeros.copy() for layer in adata.layers}


def scvi_get_latent_fields(
    mode: LatentDataType,
    scvi_latent_qzm_key: str,
    scvi_latent_qzv_key: str,
):
    """TODO add docstring"""
    if mode == "dist":
        latent_fields = [
            ObsmField(
                REGISTRY_KEYS.LATENT_QZM_KEY,
                scvi_latent_qzm_key,
            ),
            ObsmField(
                REGISTRY_KEYS.LATENT_QZV_KEY,
                scvi_latent_qzv_key,
            ),
        ]
    else:
        raise ValueError(f"Unknown latent mode: {mode}")
    latent_fields.append(
        StringUnsField(
            REGISTRY_KEYS.LATENT_MODE_KEY,
            _ADATA_LATENT_UNS_KEY,
        ),
    )
    return latent_fields
