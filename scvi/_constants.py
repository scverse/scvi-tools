from typing import NamedTuple


class _CONSTANTS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch_indices"
    LABELS_KEY: str = "labels"
    PROTEIN_EXP_KEY: str = "protein_expression"
    CAT_COVS_KEY: str = "cat_covs"
    CONT_COVS_KEY: str = "cont_covs"


_CONSTANTS = _CONSTANTS_NT()
