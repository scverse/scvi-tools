from __future__ import annotations

from functools import wraps
from typing import callable


def unsupported_if_adata_minified(fn: callable) -> callable:
    """Decorator to raise an error if the model's `adata` is minified."""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "minified_data_type", None) is not None:
            raise ValueError(
                f"The {fn.__qualname__} function currently does not support minified data."
            )
        return fn(self, *args, **kwargs)

    return wrapper
