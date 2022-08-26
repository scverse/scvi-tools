import warnings
from functools import wraps
from typing import Callable


def experimental(fn: Callable) -> Callable:
    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        warnings.warn("This is an experimental. Use with caution.")
        return fn(self, *args, **kwargs)

    return wrapper
