from typing import NamedTuple


class _MODULE_KEYS(NamedTuple):
    X_KEY: str = "x"
    # inference
    Z_KEY: str = "z"
    QZ_KEY: str = "qz"
    QZM_KEY: str = "qzm"
    QZV_KEY: str = "qzv"
    LIBRARY_KEY: str = "library"
    QL_KEY: str = "ql"
    BATCH_INDEX_KEY: str = "batch_index"
    Y_KEY: str = "y"
    CONT_COVS_KEY: str = "cont_covs"
    CAT_COVS_KEY: str = "cat_covs"
    SIZE_FACTOR_KEY: str = "size_factor"
    # generative
    PX_KEY: str = "px"
    PL_KEY: str = "pl"
    PZ_KEY: str = "pz"
    # loss
    KL_L_KEY: str = "kl_divergence_l"
    KL_Z_KEY: str = "kl_divergence_z"


MODULE_KEYS = _MODULE_KEYS()
