from typing import NamedTuple


class _DRVI_MODULE_KEYS(NamedTuple):
    # per-split decoder scale logits before aggregation, consumed by the interpretability mixin.
    # Value is a tensor of shape ``(*, n_split, n_genes)`` in inspect mode, else ``None``.
    PX_UNAGGREGATED_PARAMS_KEY: str = "px_unaggregated_params"


DRVI_MODULE_KEYS = _DRVI_MODULE_KEYS()
