"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
import warnings

from ._constants import METRIC_KEYS, REGISTRY_KEYS
from ._settings import settings

# this import needs to come after prior imports to prevent circular import
from . import data, model, external, utils, criticism

from importlib.metadata import version

package_name = "scvi-tools"
__version__ = version(package_name)

settings.verbosity = logging.INFO

# Jax sets the root logger, this prevents double output.
scvi_logger = logging.getLogger("scvi")
scvi_logger.propagate = False


__all__ = [
    "settings",
    "REGISTRY_KEYS",
    "METRIC_KEYS",
    "data",
    "model",
    "external",
    "utils",
    "criticism",
]
