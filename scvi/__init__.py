__author__ = "Romain Lopez"
__email__ = "romain_lopez@berkeley.edu"
__version__ = "0.6.7"

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._constants import _CONSTANTS
from ._settings import settings
from . import dataset, models

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())

# this prevents double outputs
logger.propagate = False

test_var = "test"

__all__ = ["settings", "_CONSTANTS", "dataset", "models"]
