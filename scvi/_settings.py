import logging


def set_global_verbosity(level: str):
    logging.getLogger("scvi").setLevel(level)
