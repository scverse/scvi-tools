import importlib
from functools import wraps
from typing import Callable


class classproperty:
    """Read-only class property decorator.

    Source: https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


def dependencies(*packages) -> Callable:
    """Checks for dependencies prior to executing a function.

    Used as a decorator, raises a `ModuleNotFoundError` if the specified
    package(s) are not installed.

    Parameters
    ----------
    packages
        A string or list of strings of packages to check for.
    """

    def decorator(fn: Callable) -> Callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            try:
                for package in packages:
                    importlib.import_module(package)
            except (ImportError, ModuleNotFoundError) as err:
                raise ImportError(
                    f"Please install {packages} to use this functionality."
                ) from err

            return fn(*args, **kwargs)

        return wrapper

    return decorator
