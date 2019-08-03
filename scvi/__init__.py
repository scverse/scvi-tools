# -*- coding: utf-8 -*-

"""Top-level package for scVI-dev."""

__author__ = "Romain Lopez"
__email__ = "romain_lopez@berkeley.edu"
__version__ = "0.4.1"

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._settings import set_verbosity

logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())

# default to INFO level logging for the scvi package
set_verbosity(logging.INFO)

__all__ = ["set_verbosity"]
