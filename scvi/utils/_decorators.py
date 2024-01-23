from functools import wraps
from typing import Callable


def unsupported_if_adata_minified(fn: Callable) -> Callable:
    """Decorator to raise an error if the model's `adata` is minified."""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "minified_data_type", None) is not None:
            raise ValueError(
                f"The {fn.__qualname__} function currently does not support minified data."
            )
        return fn(self, *args, **kwargs)

    return wrapper


class classproperty:
    """Read-only class property decorator.

    Source: https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)
