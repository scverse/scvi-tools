import logging
from typing import Optional

from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS
from scvi._types import LatentDataType
from scvi.data._constants import _ADATA_LATENT_UNS_KEY
from scvi.data.fields import ObsmField, StringUnsField

logger = logging.getLogger(__name__)

_SCVI_LATENT_SAMPLES = "_scvi_latent_samples"
_SCVI_LATENT_QZM = "_scvi_latent_qzm"
_SCVI_LATENT_QZV = "_scvi_latent_qzv"


class LatentModeMixin:
    """Functionality relevant to models that support latent mode."""

    @property
    def latent_data_type(self) -> Optional[LatentDataType]:
        return (
            self.adata_manager.get_from_registry(REGISTRY_KEYS.LATENT_MODE_KEY)
            if REGISTRY_KEYS.LATENT_MODE_KEY in self.adata_manager.data_registry
            else None
        )

    def _get_latent_adata_from_adata(
        self,
        mode: LatentDataType,
        use_latent_key: str = "X_latent",
        use_latent_qzm_key: str = "X_latent_qzm",
        use_latent_qzv_key: str = "X_latent_qzv",
    ):
        if mode == "dist":
            self.adata.obsm[_SCVI_LATENT_QZM] = self.adata.obsm[use_latent_qzm_key]
            self.adata.obsm[_SCVI_LATENT_QZV] = self.adata.obsm[use_latent_qzv_key]
        else:
            raise ValueError(f"Unknown latent mode: {mode}")
        self.adata.uns[_ADATA_LATENT_UNS_KEY] = mode
        del self.adata.raw
        all_zeros = csr_matrix(self.adata.X.shape)
        self.adata.X = all_zeros.copy()
        self.adata.layers = {layer: all_zeros.copy() for layer in self.adata.layers}

    @staticmethod
    def _get_latent_fields(mode: LatentDataType):
        if mode == "dist":
            latent_fields = [
                ObsmField(
                    REGISTRY_KEYS.LATENT_QZM_KEY,
                    _SCVI_LATENT_QZM,
                ),
                ObsmField(
                    REGISTRY_KEYS.LATENT_QZV_KEY,
                    _SCVI_LATENT_QZV,
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

    def to_latent_mode(
        self,
        mode: LatentDataType,
        use_latent_key: str = "X_latent",
        use_latent_qzm_key: str = "X_latent_qzm",
        use_latent_qzv_key: str = "X_latent_qzv",
    ):
        self._get_latent_adata_from_adata(
            mode, use_latent_key, use_latent_qzm_key, use_latent_qzv_key
        )
        self.adata_manager.register_new_fields(self.__class__._get_latent_fields(mode))
        self.module.latent_data_type = mode
