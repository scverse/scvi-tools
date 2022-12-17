import sys
from typing import Iterable, Literal

from rich.console import Console
from rich.progress import track as track_base
from tqdm import tqdm as tqdm_base

from scvi import settings


def track(
    sequence: Iterable,
    description: str = "Working...",
    disable: bool = False,
    style: Literal["rich", "tqdm"] = None,
    **kwargs,
):
    """
    Progress bar with `'rich'` and `'tqdm'` styles.

    Parameters
    ----------
    sequence
        Iterable sequence.
    description
        First text shown to left of progress bar.
    disable
        Switch to turn off progress bar.
    style
        One of ["rich", "tqdm"]. "rich" is interactive
        and is not persistent after close.
    **kwargs
        Keyword args to tqdm or rich.

    Examples
    --------
    >>> from scvi.utils import track
    >>> my_list = [1, 2, 3]
    >>> for i in track(my_list): print(i)
    """
    if style is None:
        style = settings.progress_bar_style
    if style not in ["rich", "tqdm"]:
        raise ValueError("style must be one of ['rich', 'tqdm']")
    if disable:
        return sequence
    if style == "tqdm":
        # fixes repeated pbar in jupyter
        # see https://github.com/tqdm/tqdm/issues/375
        if hasattr(tqdm_base, "_instances"):
            for instance in list(tqdm_base._instances):
                tqdm_base._decr_instances(instance)
        return tqdm_base(sequence, desc=description, file=sys.stdout, **kwargs)
    else:
        in_colab = "google.colab" in sys.modules
        force_jupyter = None if not in_colab else True
        console = Console(force_jupyter=force_jupyter)
        return track_base(sequence, description=description, console=console, **kwargs)
