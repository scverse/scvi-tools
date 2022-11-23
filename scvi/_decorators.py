from functools import wraps
from typing import Callable, List, Union


class classproperty:
    """
    Read-only class property decorator.

    Source: https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


def dependencies(packages: Union[str, List[str]]) -> Callable:
    """
    Decorator to check for dependencies.

    Parameters
    ----------
    packages
        A string or list of strings of packages to check for.
    """
    if isinstance(packages, str):
        packages = [packages]

    def decorator(fn: Callable) -> Callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            try:
                import importlib

                for package in packages:
                    importlib.import_module(package)
            except ImportError:
                raise ImportError(
                    f"Please install {packages} to use this functionality."
                )
            return fn(*args, **kwargs)

        return wrapper

    return decorator
