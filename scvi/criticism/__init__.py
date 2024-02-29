from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("sparse", "xarray")


from ._ppc import PosteriorPredictiveCheck  # noqa

__all__ = ["PosteriorPredictiveCheck"]
