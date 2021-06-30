from ._cell2location_v1 import Cell2location
from ._reference_module import RegressionModel, compute_cluster_averages

__all__ = ["Cell2location", "RegressionModel", "compute_cluster_averages"]

from pyro.distributions import constraints
from pyro.distributions.transforms import SoftplusTransform
from torch.distributions import biject_to, transform_to


@biject_to.register(constraints.positive)
@transform_to.register(constraints.positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()
