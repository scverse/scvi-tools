from collections.abc import Sequence

import numpy as np

# import xarray as xr
import pandas as pd
import torch
import torch.distributions as dist
from anndata import AnnData
from torch import Tensor
from tqdm import tqdm


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    # below params let you pass in the latent reps already, if we have already computed them
    # locs and scales must already be for the cells of the desired sample, as indices and sample
    # are ignored if locs and scales are passed in
    locs: np.ndarray | None = None,  # TODO add explanation of these vars
    scales: np.ndarray | None = None,
    sample: str | int | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
) -> dist.Distribution:
    self._check_if_trained(warn=False)

    if locs is not None and scales is not None:
        qu_loc = torch.from_numpy(locs).T
        qu_scale = torch.from_numpy(scales).T
    else:
        # TODO: If latent reps aren't passed in, I think I need to modify for if model is MrVI
        # since it's jax not pytorch, would use get_jit_inference_fn not inference
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

            outputs = self.module.inference(self.module._get_inference_input(tensors))

            qu_locs.append(outputs["qu"].loc)
            qu_scales.append(outputs["qu"].scale)

        # transpose because we need num cells to be rightmost dimension for mixture
        qu_loc = torch.cat(qu_locs, 0).T
        qu_scale = torch.cat(qu_scales, 0).T

    return dist.MixtureSameFamily(
        dist.Categorical(torch.ones(qu_loc.shape[1])), dist.Normal(qu_loc, qu_scale)
    )


# TODO: add function headers, descriptions of new params
def differential_abundance(
    self,
    adata: AnnData | None = None,
    # Below allows user to pass in latent representations directly instead
    # of using adata. Useful for simulation
    # TODO: change datatypes of below params?
    locs: np.ndarray | None = None,
    scales: np.ndarray | None = None,
    sample_id: np.ndarray | None = None,  # sample ids of each latent rep above.
    sample_cov_keys: list[str] | None = None,
    sample_subset: list[str] | None = None,
    compute_log_enrichment: bool = False,
    batch_size: int = 128,
) -> pd.DataFrame:
    adata = self._validate_anndata(adata)

    if locs is not None and scales is not None:  # if user passes in latent reps directly
        us = locs
        variances = scales
        unique_samples = np.unique(sample_id)
    else:
        # return dist so that we can also get the vars, and don't have redundantly get the latent
        # reps again in get_aggregated_posterior
        us, variances = self.get_latent_representation(
            adata, use_mean=True, give_z=False, batch_size=batch_size, return_dist=True
        )

        unique_samples = adata.obs[self.sample_key].unique()

    log_probs = []
    for sample_name in tqdm(unique_samples):
        # ap = self.get_aggregated_posterior(
        # adata=adata, sample=sample_name, batch_size=batch_size)

        # below code to prevent getting latent reps twice.
        if locs is not None and scales is not None:
            indices = np.where(sample_id == sample_name)
        else:
            indices = np.where(adata.obs[self.sample_key] == sample_name)[0]

        locs_per_sample = us[indices]
        scales_per_sample = variances[indices]
        ap = get_aggregated_posterior(self, locs=locs_per_sample, scales=scales_per_sample)

        log_probs_ = []
        n_splits = max(adata.n_obs // batch_size, 1)
        for u_rep in np.array_split(us, n_splits):
            log_probs_.append(ap.log_prob(torch.tensor(u_rep)).sum(-1, keepdims=True).cpu())

        log_probs.append(np.concatenate(log_probs_, axis=0))

    log_probs = np.concatenate(log_probs, 1)

    """coords = {
        "cell_name": adata.obs_names.to_numpy(),
        "sample": unique_samples,
    }
    data_vars = {
        "log_probs": (["cell_name", "sample"], log_probs),
    }
    log_probs_arr = xr.Dataset(data_vars, coords=coords)

    if sample_cov_keys is None or len(sample_cov_keys) == 0:
        return log_probs_arr"""

    indices = adata.obs_names.to_numpy()
    columns = unique_samples

    log_probs_df = pd.DataFrame(data=log_probs, index=indices, columns=columns)

    return log_probs_df

    # warning code that was at beginning of function
    """if sample_cov_keys is not None:
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
                )"""

    # TODO make code below into a separate function, use pandas instead of xarray,
    # user should pass in the sample info and covariates
    """sample_cov_log_probs_map = {}
    sample_cov_log_enrichs_map = {}

    for sample_cov_key in sample_cov_keys:
        # need to add sample_info or do a different way
        sample_cov_unique_values = self.sample_info[sample_cov_key].unique()
        per_val_log_probs = {}
        per_val_log_enrichs = {}

        for sample_cov_value in sample_cov_unique_values:
            cov_samples = (self.sample_info[self.sample_info[sample_cov_key] == sample_cov_value])[
                self.sample_key
            ].to_numpy()
            if sample_subset is not None:
                cov_samples = np.intersect1d(cov_samples, np.array(sample_subset))
            if len(cov_samples) == 0:
                continue

            sel_log_probs = log_probs_arr.log_probs.loc[{"sample": cov_samples}]
            val_log_probs = logsumexp(sel_log_probs, axis=1) - np.log(sel_log_probs.shape[1])
            per_val_log_probs[sample_cov_value] = val_log_probs

            if compute_log_enrichment:
                rest_samples = np.setdiff1d(unique_samples, cov_samples)
                if len(rest_samples) == 0:
                    warnings.warn(
                        f"All samples have {sample_cov_key}={sample_cov_value}. Skipping log "
                        "enrichment computation.",
                        UserWarning,
                        stacklevel=2,
                    )
                    continue

                rest_log_probs = log_probs_arr.log_probs.loc[{"sample": rest_samples}]
                rest_val_log_probs = logsumexp(rest_log_probs, axis=1) - np.log(
                    rest_log_probs.shape[1]
                )
                enrichment_scores = val_log_probs - rest_val_log_probs
                per_val_log_enrichs[sample_cov_value] = enrichment_scores

        sample_cov_log_probs_map[sample_cov_key] = DataFrame.from_dict(per_val_log_probs)
        if compute_log_enrichment and len(per_val_log_enrichs) > 0:
            sample_cov_log_enrichs_map[sample_cov_key] = DataFrame.from_dict(per_val_log_enrichs)

    coords = {
        "cell_name": adata.obs_names.to_numpy(),
        "sample": unique_samples,
        **{
            sample_cov_key: sample_cov_log_probs.columns
            for sample_cov_key, sample_cov_log_probs in sample_cov_log_probs_map.items()
        },
    }

    data_vars = {
        "log_probs": (["cell_name", "sample"], log_probs),
        **{
            f"{sample_cov_key}_log_probs": (
                ["cell_name", sample_cov_key],
                sample_cov_log_probs.values,
            )
            for sample_cov_key, sample_cov_log_probs in sample_cov_log_probs_map.items()
        },
    }

    if compute_log_enrichment:
        data_vars.update(
            {
                f"{sample_key}_log_enrichs": (
                    ["cell_name", sample_key],
                    sample_log_enrichs.values,
                )
                for sample_key, sample_log_enrichs in sample_cov_log_enrichs_map.items()
            }
        )
    return xr.Dataset(data_vars, coords=coords)"""
