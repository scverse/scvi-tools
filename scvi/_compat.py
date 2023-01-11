import warnings
from typing import Literal as _Literal


def Literal(*args, **kwargs):
    """Deprecate Literal."""
    warnings.warn(
        "Please import Literal from typing. This will be removed in the next release.",
        DeprecationWarning,
    )
    return _Literal(*args, **kwargs)
