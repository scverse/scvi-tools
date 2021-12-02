from typing import NamedTuple

################################
# scVI Manager Store Constants #
################################

_SCVI_UUID_KEY = "_scvi_uuid"

#############################
# scVI Setup Dict Constants #
#############################

_SETUP_DICT_KEY = "_scvi"
_SCVI_VERSION_KEY = "scvi_version"
_DATA_REGISTRY_KEY = "data_registry"
_CATEGORICAL_MAPPINGS_KEY = "categorical_mappings"
_SUMMARY_STATS_KEY = "summary_stats"

################################
# scVI Data Registry Constants #
################################

_DR_ATTR_NAME = "attr_name"
_DR_ATTR_KEY = "attr_key"

#######################################
# scVI Categorical Mappings Constants #
#######################################

_CM_ORIGINAL_KEY = "original_key"
_CM_MAPPING_KEY = "mapping"

#######################################
# scVI Joint Obs Categorical Mappings Constants #
#######################################

_JO_CM_MAPPINGS_KEY = "mappings"
_JO_CM_KEYS_KEY = "keys"
_JO_CM_N_CATS_PER_KEY = "n_cats_per_key"


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
