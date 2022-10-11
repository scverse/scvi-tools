"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging

# this import needs to come after prior imports to prevent circular import
from . import data, external, model, utils
from ._constants import REGISTRY_KEYS
from ._settings import settings

# https://github.com/python-poetry/poetry/pull/2366#issuecomment-652418094
# https://github.com/python-poetry/poetry/issues/144#issuecomment-623927302
try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
package_name = "scvi-tools"
__version__ = importlib_metadata.version(package_name)

settings.verbosity = logging.INFO
test_var = "test"

# Jax sets the root logger, this prevents double output.
scvi_logger = logging.getLogger("scvi")
scvi_logger.propagate = False

__all__ = ["settings", "REGISTRY_KEYS", "data", "model", "external", "utils"]
