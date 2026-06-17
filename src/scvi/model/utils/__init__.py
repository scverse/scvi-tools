from ._minification import get_minified_adata_scrna, get_minified_mudata

__all__ = ["get_minified_adata_scrna", "get_minified_mudata"]


def __getattr__(name: str):
    """Lazily provide private utilities."""
    if name in {
        "_de_core_for_annbatch",
        "_peakvi_differential_accessibility_result",
        "_run_annbatch_de_core",
        "_stream_atac_raw_counts_properties_from_dataloader",
    }:
        from . import _annbatch_de_core

        return getattr(_annbatch_de_core, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
