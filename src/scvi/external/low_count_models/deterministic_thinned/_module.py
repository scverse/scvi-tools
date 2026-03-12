"""Deterministic Thinned VAE module for ablation study."""

from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.external.low_count_models.thinned._module import ThinnedVAE
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import auto_move_data

if TYPE_CHECKING:
    from torch.distributions import Distribution


class DeterministicThinnedVAE(ThinnedVAE):
    """ThinnedVAE that uses encoder mean directly instead of sampling.

    This module extends ThinnedVAE by replacing the stochastic latent sampling
    with deterministic use of the encoder mean. This is an ablation to test whether
    the reparameterization trick's stochastic noise injection during VAE training
    is important for robustness.

    The key change is in `_regular_inference`: instead of using `z` sampled via
    `rsample()` from the posterior distribution, we use `qz.loc` (the posterior mean)
    directly, passed through `z_transformation`.

    This still keeps:
    - Thinned input data augmentation (from ThinnedVAE)
    - KL regularization term in loss (still computed, still affects encoder)

    Parameters
    ----------
    n_input
        Number of input features.
    min_library_size
        Minimum target library size for thinning. Default is 10.
        Thinned library sizes are sampled log-uniformly between this
        value and the observed library size.
    **kwargs
        Additional keyword arguments passed to
        :class:`~scvi.external.low_count_models.thinned.ThinnedVAE`.

    See Also
    --------
    :class:`~scvi.external.low_count_models.thinned.ThinnedVAE`
    :class:`~scvi.module.VAE`
    """

    @auto_move_data
    def _regular_inference(
        self,
        x: torch.Tensor,
        batch_index: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, torch.Tensor | Distribution | None]:
        """Run the regular inference process with deterministic z.

        Unlike the standard VAE, this uses the encoder mean directly
        instead of sampling from the posterior distribution.
        """
        x_ = x
        if self.use_observed_lib_size:
            library = torch.log(x.sum(1)).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log1p(x_)

        if cont_covs is not None and self.encode_covariates:
            encoder_input = torch.cat((x_, cont_covs), dim=-1)
        else:
            encoder_input = x_
        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        if self.batch_representation == "embedding" and self.encode_covariates:
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            encoder_input = torch.cat([encoder_input, batch_rep], dim=-1)
            qz, _ = self.z_encoder(encoder_input, *categorical_input)
        else:
            qz, _ = self.z_encoder(encoder_input, batch_index, *categorical_input)

        # KEY CHANGE: Use mean instead of sampling
        z = self.z_encoder.z_transformation(qz.loc)

        ql = None
        if not self.use_observed_lib_size:
            if self.batch_representation == "embedding":
                ql, library_encoded = self.l_encoder(encoder_input, *categorical_input)
            else:
                ql, library_encoded = self.l_encoder(
                    encoder_input, batch_index, *categorical_input
                )
            library = library_encoded

        if n_samples > 1:
            z = z.unsqueeze(0).expand((n_samples, z.size(0), z.size(1)))
            if self.use_observed_lib_size:
                library = library.unsqueeze(0).expand(
                    (n_samples, library.size(0), library.size(1))
                )
            else:
                library = ql.sample((n_samples,))

        return {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.QL_KEY: ql,
            MODULE_KEYS.LIBRARY_KEY: library,
        }
