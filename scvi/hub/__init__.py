from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("huggingface_hub")


from ._hub_metadata import HubMetadata, HubModelCardHelper  # noqa
from ._hub_model import HubModel  # noqa

__all__ = [
    "HubModel",
    "HubMetadata",
    "HubModelCardHelper",
]
