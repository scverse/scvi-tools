from scvi.utils import error_on_missing_dependencies

# TODO: REMOVE THIS ONCE TORCHMRVI IS READY
error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro", "xarray")

from ._model import TorchMRVI  # noqa: E402
from ._module import TorchMRVAE  # noqa: E402

__all__ = ["TorchMRVI", "TorchMRVAE"]
