from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING

import torch
from lightning.pytorch import seed_everything
from rich.console import Console
from rich.logging import RichHandler

if TYPE_CHECKING:
    from typing import Literal

scvi_logger = logging.getLogger("scvi")


class ScviConfig:
    """Config manager for scvi-tools.

    Examples
    --------
    To set the seed

    >>> scvi.settings.seed = 1

    To set the batch size for functions like `SCVI.get_latent_representation`

    >>> scvi.settings.batch_size = 1024

    To set the progress bar style, choose one of "rich", "tqdm"

    >>> scvi.settings.progress_bar_style = "rich"

    To set the verbosity

    >>> import logging
    >>> scvi.settings.verbosity = logging.INFO

    To set the number of threads PyTorch will use

    >>> scvi.settings.num_threads = 2

    To prevent Jax from preallocating GPU memory on start (default)

    >>> scvi.settings.jax_preallocate_gpu_memory = False
    """

    def __init__(
        self,
        verbosity: int = logging.INFO,
        progress_bar_style: Literal["rich", "tqdm"] = "tqdm",
        batch_size: int = 128,
        seed: int | None = None,
        logging_dir: str = "./scvi_log/",
        dl_num_workers: int = 0,
        dl_persistent_workers: bool = False,
        jax_preallocate_gpu_memory: bool = False,
        warnings_stacklevel: int = 2,
    ):
        self.warnings_stacklevel = warnings_stacklevel
        self.seed = seed
        self.batch_size = batch_size
        if progress_bar_style not in ["rich", "tqdm"]:
            raise ValueError("Progress bar style must be in ['rich', 'tqdm']")
        self.progress_bar_style = progress_bar_style
        self.logging_dir = logging_dir
        self.dl_num_workers = dl_num_workers
        self.dl_persistent_workers = dl_persistent_workers
        self._num_threads = None
        self.jax_preallocate_gpu_memory = jax_preallocate_gpu_memory
        self.verbosity = verbosity

    @property
    def batch_size(self) -> int:
        """Minibatch size for loading data into the model.

        This is only used after a model is trained. Trainers have specific
        `batch_size` parameters.
        """
        return self._batch_size

    @batch_size.setter
    def batch_size(self, batch_size: int):
        """Minibatch size for loading data into the model.

        This is only used after a model is trained. Trainers have specific
        `batch_size` parameters.
        """
        self._batch_size = batch_size

    @property
    def dl_num_workers(self) -> int:
        """Number of workers for PyTorch data loaders (Default is 0)."""
        return self._dl_num_workers

    @dl_num_workers.setter
    def dl_num_workers(self, dl_num_workers: int):
        """Number of workers for PyTorch data loaders (Default is 0)."""
        self._dl_num_workers = dl_num_workers

    @property
    def dl_persistent_workers(self) -> bool:
        """Whether to use persistent_workers in PyTorch data loaders (Default is False)."""
        return self._dl_persistent_workers

    @dl_persistent_workers.setter
    def dl_persistent_workers(self, dl_persistent_workers: bool):
        """Whether to use persistent_workers in PyTorch data loaders (Default is False)."""
        self._dl_persistent_workers = dl_persistent_workers

    @property
    def logging_dir(self) -> Path:
        """Directory for training logs (default `'./scvi_log/'`)."""
        return self._logging_dir

    @logging_dir.setter
    def logging_dir(self, logging_dir: str | Path):
        self._logging_dir = Path(logging_dir).resolve()

    @property
    def num_threads(self) -> None:
        """Number of threads PyTorch will use."""
        return self._num_threads

    @num_threads.setter
    def num_threads(self, num: int):
        """Number of threads PyTorch will use."""
        self._num_threads = num
        torch.set_num_threads(num)

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
    def seed(self, seed: int | None = None):
        """Random seed for torch and numpy."""
        if seed is None:
            self._seed = None
        else:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
            # Ensure deterministic CUDA operations for Jax (see https://github.com/google/jax/issues/13672)
            if "XLA_FLAGS" not in os.environ:
                os.environ["XLA_FLAGS"] = "--xla_gpu_deterministic_ops=true"
            else:
                os.environ["XLA_FLAGS"] += " --xla_gpu_deterministic_ops=true"
            seed_everything(seed)
            self._seed = seed

    @property
    def verbosity(self) -> int:
        """Verbosity level (default `logging.INFO`)."""
        return self._verbosity

    @verbosity.setter
    def verbosity(self, level: str | int):
        """Sets logging configuration for scvi based on chosen level of verbosity.

        If "scvi" logger has no StreamHandler, add one.
        Else, set its level to `level`.

        Parameters
        ----------
        level
            Sets "scvi" logging level to `level`
        force_terminal
            Rich logging option, set to False if piping to file output.
        """
        self._verbosity = level
        scvi_logger.setLevel(level)
        if len(scvi_logger.handlers) == 0:
            console = Console(force_terminal=True)
            if console.is_jupyter is True:
                console.is_jupyter = False
            ch = RichHandler(level=level, show_path=False, console=console, show_time=False)
            formatter = logging.Formatter("%(message)s")
            ch.setFormatter(formatter)
            scvi_logger.addHandler(ch)
        else:
            scvi_logger.setLevel(level)

    @property
    def warnings_stacklevel(self) -> int:
        """Stacklevel for warnings."""
        return self._warnings_stacklevel

    @warnings_stacklevel.setter
    def warnings_stacklevel(self, stacklevel: int):
        """Stacklevel for warnings."""
        self._warnings_stacklevel = stacklevel

    def reset_logging_handler(self):
        """Resets "scvi" log handler to a basic RichHandler().

        This is useful if piping outputs to a file.
        """
        scvi_logger.removeHandler(scvi_logger.handlers[0])
        ch = RichHandler(level=self._verbosity, show_path=False, show_time=False)
        formatter = logging.Formatter("%(message)s")
        ch.setFormatter(formatter)
        scvi_logger.addHandler(ch)

    @property
    def jax_preallocate_gpu_memory(self):
        """Jax GPU memory allocation settings.

        If False, Jax will ony preallocate GPU memory it needs.
        If float in (0, 1), Jax will preallocate GPU memory to that
        fraction of the GPU memory.
        """
        return self._jax_gpu

    @jax_preallocate_gpu_memory.setter
    def jax_preallocate_gpu_memory(self, value: float | bool):
        # see https://jax.readthedocs.io/en/latest/gpu_memory_allocation.html#gpu-memory-allocation
        if value is False:
            os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        elif isinstance(value, float):
            if value >= 1 or value <= 0:
                raise ValueError("Need to use a value between 0 and 1")
            # format is ".XX"
            os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = str(value)[1:4]
        else:
            raise ValueError("value not understood, need bool or float in (0, 1)")
        self._jax_gpu = value


settings = ScviConfig()
