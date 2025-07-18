from scvi.utils import is_package_installed

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

if is_package_installed("numpyro") and is_package_installed("jax"):
    from ._jaxvae import JaxVAE

    __all__ += ["JaxVAE"]
