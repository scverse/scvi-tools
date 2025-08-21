from scvi.utils import error_on_missing_dependencies

from ._amortizedlda import AmortizedLDAPyroModule
from ._autozivae import AutoZIVAE
from ._classifier import Classifier
from ._mrdeconv import MRDeconv
from ._multivae import MULTIVAE
from ._peakvae import PEAKVAE
from ._scanvae import SCANVAE
from ._totalvae import TOTALVAE
from ._vae import LDVAE, VAE
from ._vaec import VAEC

__all__ = [
    "VAE",
    "LDVAE",
    "TOTALVAE",
    "AutoZIVAE",
    "SCANVAE",
    "Classifier",
    "PEAKVAE",
    "VAEC",
    "MRDeconv",
    "MULTIVAE",
    "AmortizedLDAPyroModule",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxVAE":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._jaxvae import JaxVAE as _JaxVAE

        return _JaxVAE
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
