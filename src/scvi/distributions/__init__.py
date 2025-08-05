from scvi.utils import is_package_installed

from ._beta_binomial import BetaBinomial
from ._negative_binomial import (
    NegativeBinomial,
    NegativeBinomialMixture,
    Poisson,
    ZeroInflatedNegativeBinomial,
)
from ._normal import Normal

__all__ = [
    "NegativeBinomial",
    "NegativeBinomialMixture",
    "ZeroInflatedNegativeBinomial",
    "Poisson",
    "BetaBinomial",
    "Normal",
]

if is_package_installed("numpyro") and is_package_installed("jax"):
    from ._negative_binomial import JaxNegativeBinomialMeanDisp

    __all__ += ["JaxNegativeBinomialMeanDisp"]
