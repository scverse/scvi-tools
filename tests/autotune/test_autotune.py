from scvi.autotune import Autotune
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_autotune():
    adata = synthetic_iid()
    tuner = Autotune(adata=adata, model=SCVI)
    tuner.tune_scvi_asha()
