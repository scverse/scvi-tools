from ._classifier import _gaussian_process_classifier, plot_DE_results
from ._dataclass import DifferentialExpressionResults
from ._de_core import _niche_de_core
from ._de_utils import adjusted_nearest_neighbors

# from ._differential import DifferentialComputation

__all__ = [
    "_niche_de_core",
    "adjusted_nearest_neighbors",
    # "DifferentialComputation",
    # "_fdr_de_prediction",
    "_gaussian_process_classifier",
    "plot_DE_results",
    "DifferentialExpressionResults",
]
