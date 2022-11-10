from typing import NamedTuple, Optional


class _METRIC_KEYS_NT(NamedTuple):
    TRAIN: str = "train"
    VALIDATION: str = "validation"
    TEST: str = "test"
    LOSS: str = "loss"
    ELBO: str = "elbo"
    REC_LOSS: str = "reconstruction_loss"
    CLASSIFICATION_LOSS: str = "classification_loss"
    KL_LOCAL: str = "kl_local"
    KL_GLOBAL: str = "kl_global"
    SUM: str = "sum"
    N_OBS: str = "n_obs"
    N_OBS_MINIBATCH: str = "n_obs_minibatch"


METRIC_KEYS = _METRIC_KEYS_NT()


def get_metric_key(
    metric: str, mode: Optional[str] = None, metric_type: Optional[str] = None
) -> str:
    """Get metric key for logging with PyTorch Lightning."""
    assert (
        mode in [METRIC_KEYS.TRAIN, METRIC_KEYS.VALIDATION, METRIC_KEYS.TEST]
        or mode is None
    )
    assert metric_type in [METRIC_KEYS.SUM] or metric_type is None
    assert metric in [
        METRIC_KEYS.LOSS,
        METRIC_KEYS.ELBO,
        METRIC_KEYS.REC_LOSS,
        METRIC_KEYS.CLASSIFICATION_LOSS,
        METRIC_KEYS.KL_LOCAL,
        METRIC_KEYS.KL_GLOBAL,
        METRIC_KEYS.N_OBS,
        METRIC_KEYS.N_OBS_MINIBATCH,
    ]
    mode = f"{mode}_" if mode else ""
    metric_type = f"_{metric_type}" if metric_type else ""
    return mode + metric + metric_type
