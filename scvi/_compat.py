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
import sys


def track(*args, **kwargs):
    in_colab = "google.colab" in sys.modules
    force_jupyter = None if not in_colab else True
    console = Console(force_jupyter=force_jupyter)

    if "disable" in kwargs.keys():
        disable = kwargs.pop("disable")
    else:
        disable = False
    if disable:
        return args[0]
    else:
        return track_base(*args, **kwargs, console=console)
