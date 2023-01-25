import scvi
from scvi.autotune._defaults import DEFAULTS


def test_default_search_spaces():
    # test that all default search spaces are valid
    for model_cls in DEFAULTS:
        search_space = DEFAULTS[model_cls]
        manager = scvi.autotune.TunerManager(model_cls)
        registry = manager._registry["tunables"]
        for key in search_space:
            assert key in registry
