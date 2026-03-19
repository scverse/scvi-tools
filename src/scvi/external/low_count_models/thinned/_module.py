"""Thinned VAE module for ablation study."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scvi import REGISTRY_KEYS
from scvi.external.low_count_models.utils import binomial_split, sample_thinning_probs
from scvi.module._constants import MODULE_KEYS
from scvi.module._vae import VAE

if TYPE_CHECKING:
    import torch


class ThinnedVAE(VAE):
    """VAE that trains on randomly thinned data.

    This module extends the standard VAE by training on thinned versions of the
    input data. During training, each batch is randomly thinned using binomial
    thinning with per-cell probabilities sampled to produce target library sizes
    that are log-uniform between min_library_size and the observed library size.

    This is an ablation model to test whether exposure to thinned data during
    training (affecting encoder weights and BatchNorm statistics) is sufficient
    for robustness to low-UMI cells, or whether explicit invariance supervision
    (like CCO loss in JointEmbeddingSCVI) is required.

    Parameters
    ----------
    n_input
        Number of input features.
    min_library_size
        Minimum target library size for thinning. Default is 10.
        Thinned library sizes are sampled log-uniformly between this
        value and the observed library size.
    **kwargs
        Additional keyword arguments passed to :class:`~scvi.module.VAE`.

    See Also
    --------
    :class:`~scvi.module.VAE`
    :class:`~scvi.external.JointEmbeddingSCVI`
    """

    def __init__(
        self,
        n_input: int,
        *args,
        min_library_size: float = 10.0,
        **kwargs,
    ):
        super().__init__(n_input, *args, **kwargs)
        self.min_library_size = min_library_size

    def _get_inference_input(
        self,
        tensors: dict[str, torch.Tensor | None],
        full_forward_pass: bool = False,
    ) -> dict[str, torch.Tensor | None]:
        """Get input tensors for the inference process.

        During training, the input data is randomly thinned before being passed
        to the encoder. At inference time, the original data is used.
        """
        if full_forward_pass or self.minified_data_type is None:
            loader = "full_data"
        elif self.minified_data_type in [
            "latent_posterior_parameters",
            "latent_posterior_parameters_with_counts",
        ]:
            loader = "minified_data"
        else:
            raise NotImplementedError(f"Unknown minified-data type: {self.minified_data_type}")

        if loader == "minified_data":
            return {
                MODULE_KEYS.QZM_KEY: tensors[REGISTRY_KEYS.LATENT_QZM_KEY],
                MODULE_KEYS.QZV_KEY: tensors[REGISTRY_KEYS.LATENT_QZV_KEY],
                REGISTRY_KEYS.OBSERVED_LIB_SIZE: tensors[REGISTRY_KEYS.OBSERVED_LIB_SIZE],
            }

        x = tensors[REGISTRY_KEYS.X_KEY]

        if self.training:
            p = sample_thinning_probs(x, min_library_size=self.min_library_size)
            x, _ = binomial_split(x, p=p)

        return {
            MODULE_KEYS.X_KEY: x,
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            MODULE_KEYS.CONT_COVS_KEY: tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None),
            MODULE_KEYS.CAT_COVS_KEY: tensors.get(REGISTRY_KEYS.CAT_COVS_KEY, None),
        }
