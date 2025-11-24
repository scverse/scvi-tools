from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")

from ._model import JaxMRVI  # noqa: E402
from ._module import JaxMRVAE  # noqa: E402

__all__ = ["JaxMRVI", "JaxMRVAE"]
