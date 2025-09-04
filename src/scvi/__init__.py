"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging

from ._constants import REGISTRY_KEYS
from ._settings import settings

# this import needs to come after prior imports to prevent circular import
from . import data, model, external, utils

from importlib.metadata import version
from scvi.utils import error_on_missing_dependencies

package_name = "scvi-tools"
__version__ = version(package_name)

settings.verbosity = logging.INFO

# Jax sets the root logger, this prevents double output.
scvi_logger = logging.getLogger("scvi")
scvi_logger.propagate = False


__all__ = [
    "settings",
    "REGISTRY_KEYS",
    "data",
    "model",
    "external",
    "utils",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "criticism":
        error_on_missing_dependencies("xarray", "sparse")
        from . import criticism as _criticism

        return _criticism
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
