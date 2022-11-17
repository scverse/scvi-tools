import multiprocessing

import pytest

import scvi


@pytest.mark.optional
def test_basic():
    curr_proc = multiprocessing.current_process()
    print("current process: ", curr_proc.name, curr_proc._identity)
    model_cls = scvi.model.SCVI
    scvi.autotune.ModelTuner(model_cls)
    print("current process: ", curr_proc.name, curr_proc._identity)
