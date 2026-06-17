from ._minification import get_minified_adata_scrna, get_minified_mudata

__all__ = ["get_minified_adata_scrna", "get_minified_mudata"]


def __getattr__(name: str):
    """Lazily provide private utilities."""
    if name == "_de_core_for_annbatch":
        from ._annbatch_de_core import _de_core_for_annbatch

        return _de_core_for_annbatch
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
