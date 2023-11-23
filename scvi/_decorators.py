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


def dependencies(*modules) -> Callable:
    """Decorator to check for dependencies.

    Parameters
    ----------
    packages
        A string or list of strings of packages to check for.
    """

    def decorator(fn: Callable) -> Callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            missing_modules = []
            for module in modules:
                try:
                    importlib.import_module(module)
                except ImportError:
                    missing_modules.append(module)

            if len(missing_modules) > 0:
                raise ModuleNotFoundError(
                    f"Please install {missing_modules} to use this functionality."
                )
            return fn(*args, **kwargs)

        return wrapper

    return decorator
