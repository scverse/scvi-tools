from typing import NamedTuple


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


SAVE_KEYS = _SAVE_KEYS_NT()
