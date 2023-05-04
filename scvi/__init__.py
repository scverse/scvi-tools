"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging

from ._constants import METRIC_KEYS, REGISTRY_KEYS
from ._settings import settings

# this import needs to come after prior imports to prevent circular import
from . import autotune, data, model, external, utils

from importlib.metadata import version

package_name = "scvi-tools"
__version__ = version(package_name)

settings.verbosity = logging.INFO
test_var = "test"

# Jax sets the root logger, this prevents double output.
scvi_logger = logging.getLogger("scvi")
scvi_logger.propagate = False

__all__ = [
    "settings",
    "REGISTRY_KEYS",
    "METRIC_KEYS",
    "autotune",
    "data",
    "model",
    "external",
    "utils",
]
