import pytest

import scvi


@pytest.mark.optional
def test_basic():
    model_cls = scvi.model.SCVI
    scvi.autotune.ModelTuner(model_cls)
