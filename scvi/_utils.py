import sys
import warnings
from textwrap import dedent
from typing import Iterable

from tqdm import tqdm as tqdm_base
from tqdm.rich import tqdm as tqdm_rich

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return tqdm_rich(sequence, desc=description, **kwargs)


def _doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj

    return dec
