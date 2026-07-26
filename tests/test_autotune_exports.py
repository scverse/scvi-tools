import pytest


def test_exports_the_real_implementations():
    """``scvi.autotune`` must re-export the implementations, never a placeholder.

    The placeholders exist only so that importing the package raises a helpful error when the
    optional dependencies are missing. When they are installed, shadowing the real classes hides
    their attributes and breaks the API documentation.
    """
    pytest.importorskip("hyperopt")
    pytest.importorskip("ray.tune")

    import scvi.autotune

    assert scvi.autotune.AutotuneExperiment.__module__ == "scvi.autotune._experiment"
    assert scvi.autotune.ScibTuneReportCheckpointCallback.__module__ == "scvi.autotune._experiment"
    assert scvi.autotune.run_autotune.__module__ == "scvi.autotune._tune"
    assert hasattr(scvi.autotune.AutotuneExperiment, "result_grid")
