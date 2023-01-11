import warnings
from typing import Literal as _Literal


class Literal(_Literal):
    """Shim Literal."""

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "Please import Literal from typing. This will be removed in the next release.",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)
