"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._constants import _CONSTANTS
from ._settings import settings

# this import needs to come after prior imports to prevent circular import
from . import data, model

# https://github.com/python-poetry/poetry/pull/2366#issuecomment-652418094
# https://github.com/python-poetry/poetry/issues/144#issuecomment-623927302
try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
package_name = "scvi-tools"
__version__ = importlib_metadata.version(package_name)

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())

# this prevents double outputs
logger.propagate = False

test_var = "test"

__all__ = ["settings", "_CONSTANTS", "data", "model"]
