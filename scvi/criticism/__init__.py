from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("xarray", "sparse")


from ._ppc import PosteriorPredictiveCheck  # noqa

__all__ = ["PosteriorPredictiveCheck"]
