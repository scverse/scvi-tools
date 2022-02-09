from typing import NamedTuple

# scVI Manager Store Constants
# ----------------------------
# Keys for UUIDs used for referencing model class manager stores.

_SCVI_UUID_KEY = "_scvi_uuid"
_MANAGER_UUID_KEY = "_scvi_manager_uuid"

# scVI Registry Constants
# -----------------------
# Keys used in the scVI registry.

_SCVI_VERSION_KEY = "scvi_version"
_MODEL_NAME_KEY = "model_name"
_SETUP_KWARGS_KEY = "setup_kwargs"
_FIELD_REGISTRIES_KEY = "field_registries"
_DATA_REGISTRY_KEY = "data_registry"
_STATE_REGISTRY_KEY = "state_registry"
_SUMMARY_STATS_KEY = "summary_stats"

# scVI Data Registry Constants
# ----------------------------
# Keys used in the data registry.

_DR_ATTR_NAME = "attr_name"
_DR_ATTR_KEY = "attr_key"

# AnnData Object Constants
# ------------------------
# AnnData object attribute names.


class _ADATA_ATTRS_NT(NamedTuple):
    X: str = "X"
    LAYERS: str = "layers"
    OBS: str = "obs"
    OBSM: str = "obsm"
    VAR: str = "var"
    VARM: str = "varm"


_ADATA_ATTRS = _ADATA_ATTRS_NT()
