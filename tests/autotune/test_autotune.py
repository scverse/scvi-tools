from ray.tune import loguniform

from scvi.autotune import Autotune
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_autotune():
    adata = synthetic_iid()
    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
        "kl_global_validation",
    ]
    config = {"dropout_rate": loguniform(1e-4, 1e-1)}
    tuner = Autotune(adata, SCVI, "scvi_tuner", "elbo_validation", metrics, config)
    tuner()
