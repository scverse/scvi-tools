import logging
from typing import Union

scvi_logger = logging.getLogger("scvi")


def set_verbosity(level: Union[str, int]):
    """Sets logging configuration for scvi based on chosen level of verbosity."""
    scvi_logger.setLevel(level)
    has_streamhandler = False
    for handler in scvi_logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            has_streamhandler = True
    if not has_streamhandler:
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s - %(name)s | %(message)s"
        )
        ch.setFormatter(formatter)
        scvi_logger.addHandler(ch)
