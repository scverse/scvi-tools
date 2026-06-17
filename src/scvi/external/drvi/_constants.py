from typing import NamedTuple


class _DRVI_MODULE_KEYS(NamedTuple):
    # per-split decoder parameters before aggregation, exposed in inspect mode and consumed by the
    # interpretability mixin. Value is a dict, e.g. ``{"mean": tensor}`` of shape
    # ``(n_obs, n_split, n_genes)``.
    PX_UNAGGREGATED_PARAMS_KEY: str = "px_unaggregated_params"


DRVI_MODULE_KEYS = _DRVI_MODULE_KEYS()
