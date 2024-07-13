import warnings
from collections.abc import Sequence

import numpy as np
import torch
import xarray as xr
from anndata import AnnData
from torch import Tensor
from torch.distributions import Distribution, Normal


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    sample: str | int | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
) -> Distribution:
    self._check_if_trained(warn=False)
    adata = self._validate_anndata(adata)

    if indices is None:
        indices = np.arrange(self.adata.n_obs)
    if sample is not None:
        indices = np.intersect1d(
            np.array(indices), np.where(adata.obs[self.sample_key] == sample)[0]
        )

    dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

    qu_locs: list[Tensor] = []
    qu_scales: list[Tensor] = []
    for tensors in dataloader:
        """ Not sure how to get the u estimated posteriors.
        I don't think the following would work because the location of the
        code would be in _vaemixin.py, and I think this would work if it
        was still in mrvi/_model.py"""

        """Also I think this would return a jax distribution.
        Do I need to create a new inference fn?"""
        outputs = self.module.inference(self.module._get_inference_input(tensors))

        qu_locs.append(outputs["qu"].loc)
        qu_scales.append(outputs["qu"].scale)

    # transpose because we need num cells to be rightmost dimension for mixture
    qu_loc = torch.cat(qu_locs, 0).T
    qu_scale = torch.cat(qu_scales, 0).T

    return Distribution.MixtureSameFamily(
        Distribution.Categorical(torch.ones(qu_loc.shape[1])), Normal(qu_loc, qu_scale)
    )


def differential_abundance(
    self,
    adata: AnnData | None = None,
    sample_cov_keys: list[str] | None = None,
    batch_size: int = 128,
) -> xr.Dataset:
    adata = self._validate_anndata(adata)

    if sample_cov_keys is not None:
        for key in sample_cov_keys:
            n_cov_values = len(adata.obs[key].unique())
            n_samples = len(adata.obs[self.sample_key].unique())
            if n_cov_values > n_samples / 2:
                warnings.warn(
                    f"The covariate '{key}' does not seem to refer to a discrete key. "
                    f"It has {len(n_cov_values)} unique values, which exceeds one half of the "
                    f"total samples ({n_samples}).",
                    UserWarning,
                    stacklevel=2,
                )

    """Same issue as with get_aggregated_posterior. Not sure how I should get the u
    latent representation as get_latent_representation in _vaemixin only has the z
    representation, and get_latent_representation for mrvi uses jax."""
    # us = self.get_latent_representation(
    # adata, use_mean=True, give_z=False, batch_size=batch_size
    # )

    # log_probs = []
    # unique_samples = adata.obs[self.sample_key].unique()
    # for sample_name in tqdm(unique_samples):
    # ap = self.get_aggregated_posterior()
