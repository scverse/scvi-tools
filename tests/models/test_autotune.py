from scvi.train import Autotune
from scvi.model import SCVI
from scvi.data import synthetic_iid

def test_autotune():
    adata = synthetic_iid()
    tuner = Autotune(adata=adata, model=SCVI,n_samples = 10)
    tuner.tune_scvi_asha()