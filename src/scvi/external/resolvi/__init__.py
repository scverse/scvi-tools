from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("pyro")

from ._model import RESOLVI  # noqa: E402
from ._module import RESOLVAE  # noqa: E402

__all__ = ["RESOLVAE", "RESOLVI"]
