try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


from rich.progress import track as track_base
from rich.console import Console
from tqdm import tqdm as tqdm_base
from typing import Iterable
import sys


def track(
    sequence: Iterable,
    description: str = "Working...",
    disable: bool = False,
    style: Literal["rich", "tqdm"] = "tqdm",
    **kwargs
):
    """Progress bar with `'rich'` and `'tqdm'` styles."""

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
