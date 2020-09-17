"""scvi-tools."""

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._constants import _CONSTANTS
from ._settings import settings
from . import data, model

import toml
from pathlib import Path


def get_version(source_file):
    # https://github.com/rominf/poetry-version/blob/master/poetry_version/__init__.py
    d = Path(source_file)
    result = None
    while d.parent != d and result is None:
        d = d.parent
        pyproject_toml_path = d / "pyproject.toml"
        if pyproject_toml_path.exists():
            with open(file=str(pyproject_toml_path)) as f:
                pyproject_toml = toml.loads(f.read())
                if "tool" in pyproject_toml and "poetry" in pyproject_toml["tool"]:
                    # noinspection PyUnresolvedReferences
                    result = pyproject_toml["tool"]["poetry"]["version"]
    return result


__version__ = get_version(__file__)

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())

# this prevents double outputs
logger.propagate = False

test_var = "test"

__all__ = ["settings", "_CONSTANTS", "data", "model"]
