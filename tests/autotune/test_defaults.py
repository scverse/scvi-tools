def test_default_search_spaces():
    from scvi.autotune import TunerManager
    from scvi.autotune._defaults import DEFAULTS

    # test that all default search spaces are valid
    for model_cls in DEFAULTS:
        search_space = DEFAULTS[model_cls]
        manager = TunerManager(model_cls)
        registry = manager._registry["tunables"]
        for key in search_space:
            assert key in registry
