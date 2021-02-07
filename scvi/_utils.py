import sys
from textwrap import dedent
from typing import Iterable

from rich.console import Console
from rich.progress import track as track_base
from tqdm import tqdm as tqdm_base

from scvi import settings
from scvi._compat import Literal


def track(
    sequence: Iterable,
    description: str = "Working...",
    disable: bool = False,
    style: Literal["rich", "tqdm"] = None,
    **kwargs
):
    """Progress bar with `'rich'` and `'tqdm'` styles."""
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


def _doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj

    return dec
