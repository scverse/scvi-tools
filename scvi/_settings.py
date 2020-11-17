import logging
from typing import Union

import numpy as np
import torch
from rich.console import Console
from rich.logging import RichHandler

from ._compat import Literal

logger = logging.getLogger(__name__)
scvi_logger = logging.getLogger("scvi")


class ScviConfig:
    """
    Config manager for scvi-tools.

    Examples
    --------
    >>> import scvi
    >>> scvi.settings.seed = 1
    """

    def __init__(
        self,
        verbosity: int = logging.INFO,
        progress_bar_style: Literal["rich", "tqdm"] = "tqdm",
        batch_size: int = 128,
        seed: int = 0,
    ):

        self.verbosity = verbosity
        self.seed = seed
        self.batch_size = batch_size
        if progress_bar_style not in ["rich", "tqdm"]:
            raise ValueError("Progress bar style must be in ['rich', 'tqdm']")
        self.progress_bar_style = progress_bar_style

    @property
    def batch_size(self) -> int:
        """
        Minibatch size for loading data into the model.

        This is only used after a model is trained. Trainers have specific
        `batch_size` parameters.
        """
        return self._batch_size

    @batch_size.setter
    def batch_size(self, batch_size: int):
        """
        Minibatch size for loading data into the model.

        This is only used after a model is trained. Trainers have specific
        `batch_size` parameters.
        """
        self._batch_size = batch_size

    @property
    def progress_bar_style(self) -> str:
        """Library to use for progress bar."""
        return self._pbar_style

    @progress_bar_style.setter
    def progress_bar_style(self, pbar_style: Literal["tqdm", "rich"]):
        """Library to use for progress bar."""
        self._pbar_style = pbar_style

    @property
    def seed(self) -> int:
        """Random seed for torch and numpy."""
        return self._seed

    @seed.setter
    def seed(self, seed: int):
        """Random seed for torch and numpy."""
        torch.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        np.random.seed(seed)
        self._seed = seed

    @property
    def verbosity(self) -> int:
        """Verbosity level (default `logging.INFO`)."""
        return self._verbosity

    @verbosity.setter
    def verbosity(self, level: Union[str, int]):
        """
        Sets logging configuration for scvi based on chosen level of verbosity.

        Sets "scvi" logging level to `level`
        If "scvi" logger has no StreamHandler, add one.
        Else, set its level to `level`.
        """
        self._verbosity = level
        scvi_logger.setLevel(level)
        has_streamhandler = False
        for handler in scvi_logger.handlers:
            if isinstance(handler, RichHandler):
                handler.setLevel(level)
                logger.info(
                    "'scvi' logger already has a StreamHandler, set its level to {}.".format(
                        level
                    )
                )
                has_streamhandler = True
        if not has_streamhandler:
            console = Console(force_terminal=True)
            if console.is_jupyter is True:
                console.is_jupyter = False
            ch = RichHandler(show_path=False, console=console, show_time=False)
            formatter = logging.Formatter("%(message)s")
            ch.setFormatter(formatter)
            scvi_logger.addHandler(ch)
            logger.debug("Added StreamHandler with custom formatter to 'scvi' logger.")


settings = ScviConfig()
