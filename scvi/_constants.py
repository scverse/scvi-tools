from typing import NamedTuple


class _CONSTANTS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    LABELS_KEY: str = "labels"
    PROTEIN_EXP_KEY: str = "protein_expression"
    CAT_COVS_KEY: str = "extra_categoricals"
    CONT_COVS_KEY: str = "extra_continuous"


_CONSTANTS = _CONSTANTS_NT()
