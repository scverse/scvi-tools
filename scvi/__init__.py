"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._constants import _CONSTANTS
from ._settings import settings
from . import dataset, models

import toml
from pathlib import Path


def get_version():
    path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = toml.loads(open(str(path)).read())
    return pyproject["tool"]["poetry"]["version"]


__version__ = get_version()

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())

# this prevents double outputs
logger.propagate = False

test_var = "test"

__all__ = ["settings", "_CONSTANTS", "dataset", "models"]
