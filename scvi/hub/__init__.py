from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("huggingface_hub")


from ._metadata import HubMetadata, HubModelCardHelper  # noqa
from ._model import HubModel  # noqa

__all__ = [
    "HubModel",
    "HubMetadata",
    "HubModelCardHelper",
]
