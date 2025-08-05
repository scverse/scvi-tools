from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("pyro")

from ._model import Decipher  # noqa: E402
from ._module import DecipherPyroModule  # noqa: E402

__all__ = ["Decipher", "DecipherPyroModule"]
