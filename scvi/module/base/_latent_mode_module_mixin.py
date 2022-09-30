import logging
from typing import Optional

from scvi._types import LatentDataType
from scvi.module.base import auto_move_data

logger = logging.getLogger(__name__)


class LatentModeModuleMixin:
    """Functionality relevant to modules that support latent mode."""

    @property
    def latent_data_type(self) -> Optional[LatentDataType]:
        return self._latent_data_type

    @latent_data_type.setter
    def latent_data_type(self, latent_data_type):
        self._latent_data_type = latent_data_type

    @auto_move_data
    def inference(self, *args, **kwargs):
        if self.latent_data_type is None:
            return self._regular_inference(*args, **kwargs)
        else:
            return self._inference_no_encode(*args, **kwargs)
