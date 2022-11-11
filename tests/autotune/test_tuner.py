import scvi
from scvi import autotune


def test_basic():
    model_cls = scvi.model.SCVI
    autotune.ModelTuner(model_cls)
