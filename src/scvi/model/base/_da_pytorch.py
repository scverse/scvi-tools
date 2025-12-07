from collections.abc import Sequence

import numpy as np
import pandas as pd
import torch
import torch.distributions as dist
from anndata import AnnData
from tqdm import tqdm


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    sample: str | int | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
    dof: float | None = 3.0,
) -> dist.Distribution:
    """Compute the aggregated posterior over the ``u`` latent representations.

    Parameters
    ----------
    adata
        AnnData object to use. Defaults to the AnnData object used to initialize the model.
    sample
        Name or index of the sample to filter on. If ``None``, uses all cells.
    indices
        Indices of cells to use.
    batch_size
        Batch size to use for computing the latent representation.
    dof
        Degrees of freedom for the Student's t-distribution components.
        If ``None``, components are Normal.

    Returns
    -------
    A mixture distribution of the aggregated posterior.
    """
    self._check_if_trained(warn=False)
    adata = self._validate_anndata(adata)

    if indices is None:
        indices = np.arange(adata.n_obs)
    if sample is not None:
        indices = np.intersect1d(
            np.array(indices), np.where(adata.obs[self.sample_key] == sample)[0]
        )

    dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
    qu_loc, qu_scale = self.get_latent_representation(
        batch_size=batch_size, return_dist=True, dataloader=dataloader, give_mean=True
    )

    qu_loc = torch.tensor(qu_loc, device=self.device)  # (n_cells, n_latent_u)
    qu_scale = torch.tensor(qu_scale, device=self.device)

    if dof is None:
        components = dist.Normal(qu_loc, qu_scale)
    else:
        components = dist.StudentT(dof, qu_loc, qu_scale)
    return dist.MixtureSameFamily(
        dist.Categorical(logits=torch.ones(qu_loc.shape[0], device=self.device)),
        dist.Independent(components, 1),
    )


def differential_abundance(
    self,
    adata: AnnData | None = None,
    sample_key: str | None = None,
    batch_size: int = 128,
    num_cells_posterior: int | None = None,
    dof: float | None = None,
):
    """Compute the differential abundance between samples.

    Computes the log probabilities of each sample conditioned on the estimated
    aggregate posterior distribution of each cell.

    Parameters
    ----------
    adata
        The data object to compute the differential abundance for.
        For very large datasets, this should be a subset of the original data object.
    sample_key
        Key for the sample covariate.
    batch_size
        Minibatch size for computing the differential abundance.
    num_cells_posterior
        Maximum number of cells used to compute aggregated posterior for each sample.
    dof
        Degrees of freedom for the Student's t-distribution components for aggregated posterior.
        If ``None``, components are Normal.
    """
    adata = self._validate_anndata(adata)

    # In case user passes in a subset of model's anndata
    adata_dataloader = self._make_data_loader(adata=adata, batch_size=batch_size)
    us = self.get_latent_representation(
        batch_size=batch_size, dataloader=adata_dataloader, give_mean=True
    )
    dataloader = torch.utils.data.DataLoader(us, batch_size=batch_size)
    unique_samples = adata.obs[sample_key].unique()

    log_probs = []
    for sample_name in tqdm(unique_samples):
        indices = np.where(adata.obs[sample_key] == sample_name)[0]
        if num_cells_posterior is not None and num_cells_posterior < indices.shape[0]:
            indices = np.random.choice(indices, num_cells_posterior, replace=False)

        ap = get_aggregated_posterior(
            self, adata=adata, indices=indices, dof=dof, batch_size=batch_size
        )
        log_probs_ = []
        for u_rep in dataloader:
            u_rep = u_rep.to(self.device)
            log_probs_.append(ap.log_prob(u_rep))
        log_probs.append(torch.cat(log_probs_, axis=0).cpu().numpy())

    log_probs = np.array(log_probs).T
    log_probs_df = pd.DataFrame(data=log_probs, index=adata.obs_names, columns=unique_samples)

    adata.obsm["da_log_probs"] = log_probs_df
