from .classifier import Classifier
from .scanvi import SCANVI
from .vae import VAE, LDVAE
from .autozivae import AutoZIVAE
from .vaec import VAEC
from .jvae import JVAE
from .totalvi import TOTALVI


__doc__ = """\
scvi.models
===========

The central class:

.. autosummary::
   :toctree: .

   VAE
   TOTALVI

"""

__all__ = [
    "SCANVI",
    "VAEC",
    "VAE",
    "LDVAE",
    "JVAE",
    "Classifier",
    "AutoZIVAE",
    "TOTALVI",
]
