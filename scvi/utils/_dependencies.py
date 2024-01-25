from __future__ import annotations

import importlib
from functools import wraps
from typing import callable


def error_on_missing_dependencies(*modules):
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


def dependencies(*modules) -> callable:
    """Decorator to check for dependencies."""

    def decorator(fn: callable) -> callable:
        @wraps(fn)
        def wrapper(*args, **kwargs):
            error_on_missing_dependencies(*modules)
            return fn(*args, **kwargs)

        return wrapper

    return decorator
