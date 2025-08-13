from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("xarray", "sparse")

from ._create_criticism_report import create_criticism_report  # noqa: E402
from ._ppc import PosteriorPredictiveCheck  # noqa: E402

__all__ = ["PosteriorPredictiveCheck", "create_criticism_report"]
