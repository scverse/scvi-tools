from typing import NamedTuple


class _METHYLVI_REGISTRY_KEYS_NT(NamedTuple):
    MC_KEY: str = "mc"
    COV_KEY: str = "cov"
    ANNDATA_MODALITY_PLACEHOLDER: str = "methylation"


METHYLVI_REGISTRY_KEYS = _METHYLVI_REGISTRY_KEYS_NT()
