from typing import NamedTuple


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    LABELS_KEY: str = "labels"
    PROTEIN_EXP_KEY: str = "proteins"
    CAT_COVS_KEY: str = "extra_categorical_covs"
    CONT_COVS_KEY: str = "extra_continuous_covs"
    INDICES_KEY: str = "ind_x"
    SIZE_FACTOR_KEY: str = "size_factor"
    MINIFY_TYPE_KEY: str = "minify_type"
    LATENT_QZM_KEY: str = "latent_qzm"
    LATENT_QZV_KEY: str = "latent_qzv"
    OBSERVED_LIB_SIZE: str = "observed_lib_size"


class _SAVE_KEYS_NT(NamedTuple):
    ADATA_FNAME: str = "adata.h5ad"
    MDATA_FNAME: str = "mdata.h5mu"
    MODEL_FNAME: str = "model.pt"
    MODEL_STATE_DICT_KEY: str = "model_state_dict"
    VAR_NAMES_KEY: str = "var_names"
    ATTR_DICT_KEY: str = "attr_dict"
    LEGACY_MODEL_FNAME: str = "model_params.pt"
    LEGACY_VAR_NAMES_FNAME: str = "var_names.csv"
    LEGACY_SETUP_DICT_FNAME: str = "attr.pkl"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()
SAVE_KEYS = _SAVE_KEYS_NT()
