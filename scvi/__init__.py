# -*- coding: utf-8 -*-

"""Top-level package for scVI-dev."""

__author__ = "Romain Lopez"
__email__ = "romain_lopez@berkeley.edu"
__version__ = '0.3.0'

# Set default logging handler to avoid logging with logging.lastResort logger.
import logging
from logging import NullHandler

from ._settings import set_global_verbosity

logging.getLogger(__name__).addHandler(NullHandler())

__all__ = ["set_global_verbosity"]
