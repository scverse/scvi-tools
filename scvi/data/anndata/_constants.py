from typing import NamedTuple

################################
# scVI Manager Store Constants #
################################

_SCVI_UUID_KEY = "_scvi_uuid"
_SOURCE_SCVI_UUID_KEY = "_source_scvi_uuid"

#############################
# scVI Registry Constants #
#############################

_SCVI_VERSION_KEY = "scvi_version"
_MODEL_NAME_KEY = "model_name"
_SETUP_KWARGS_KEY = "setup_kwargs"
_FIELD_REGISTRIES_KEY = "field_registries"
_DATA_REGISTRY_KEY = "data_registry"
_STATE_REGISTRY_KEY = "state_registry"
_SUMMARY_STATS_KEY = "summary_stats"

################################
# scVI Data Registry Constants #
################################

_DR_ATTR_NAME = "attr_name"
_DR_ATTR_KEY = "attr_key"

############################
# AnnData Object Constants #
############################


class _ADATA_ATTRS_NT(NamedTuple):
    X: str = "X"
    LAYERS: str = "layers"
    OBS: str = "obs"
    OBSM: str = "obsm"
    VAR: str = "var"
    VARM: str = "varm"


_ADATA_ATTRS = _ADATA_ATTRS_NT()
