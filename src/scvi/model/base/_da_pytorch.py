from collections.abc import Sequence

import numpy as np

import pandas as pd
import torch
import torch.distributions as dist
from anndata import AnnData
from torch import Tensor
from tqdm import tqdm
from scipy.sparse import coo_matrix


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    sample: str | int | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
    dof: float | None = 3.,
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
        Degrees of freedom for the Student's t-distribution components. If ``None``, components are Normal.

    Returns
    -------
    A mixture distribution of the aggregated posterior.
    """
    self._check_if_trained(warn=False)
    adata = self._validate_anndata(adata)

    if indices is None:
        indices = np.arange(self.adata.n_obs)
    if sample is not None:
        indices = np.intersect1d(
            np.array(indices), np.where(adata.obs[self.sample_key] == sample)[0]
        )

    dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
    qu_loc, qu_scale = self.get_latent_representation(batch_size=batch_size, return_dist=True, dataloader=dataloader, give_mean=True)

    qu_loc = torch.tensor(qu_loc, device=self.device) # (n_cells, n_latent_u)
    qu_scale = torch.tensor(qu_scale, device=self.device)
    
    if dof is None:
        components = dist.Normal(qu_loc, qu_scale)
    else:
        components = dist.StudentT(dof, qu_loc, qu_scale)
    return dist.MixtureSameFamily(
        dist.Categorical(logits=torch.ones(qu_loc.shape[0], device=self.device)), dist.Independent(components, 1))

def differential_abundance(
    self,
    adata: AnnData | None = None,
    sample_key: str | None = None,
    batch_size: int = 128,
    downsample: int | None = None,
    num_cells_posterior: int | None = None,
    dof: float | None = None,
) -> pd.DataFrame:
    """Compute the differential abundance between samples.

    Computes the logarithm of the ratio of the probabilities of each sample conditioned on the
    estimated aggregate posterior distribution of each cell.

    Parameters
    ----------
    adata
        The data object to compute the differential abundance for.
    sample_key
        Key for the sample covariate.
    batch_size
        Minibatch size for computing the differential abundance.
    downsample
        Number of cells for which we will compute and store log probabilities.
    num_cells_posterior
        Maximum number of cells used to compute aggregated posterior for each sample.
    dof
        Degrees of freedom for the Student's t-distribution components for aggregated posterior. If ``None``, components are Normal.

    Returns
    -------
    DataFrame of shape (n_cells, n_samples) containing the log probabilities
    for each cell across samples. The rows correspond to cell names from `adata.obs_names`,
    and the columns correspond to unique sample identifiers.
    """

    adata = self._validate_anndata(adata)

    us = self.get_latent_representation(
        batch_size=batch_size, return_dist=False, give_mean=True
    )

    # Are we computing for all cells?
    downsampling = downsample and (downsample < us.shape[0])
    if downsampling:
        downsample_indices = np.random.choice(range(us.shape[0]), downsample, replace=False)
    else:
        downsample_indices = np.array(list(range(us.shape[0])))

    dataloader = torch.utils.data.DataLoader(us[downsample_indices], batch_size=batch_size)
    unique_samples = adata.obs[sample_key].unique()

    log_probs = []
    for sample_name in tqdm(unique_samples):
        indices = np.where(adata.obs[sample_key] == sample_name)[0]
        if num_cells_posterior is not None and num_cells_posterior < indices.shape[0]:
            indices = np.random.choice(indices, num_cells_posterior, replace=False)

        ap = get_aggregated_posterior(self, adata=adata, indices=indices, dof=dof)
        log_probs_ = []
        for u_rep in dataloader:
            u_rep = u_rep.to(self.device)
            log_probs_.append(ap.log_prob(u_rep))
        log_probs.append(torch.cat(log_probs_, axis=0).cpu().numpy())

    # If not downsampling, this is the full (num_cells, num_samples) array
    # If we are downsampling this is (downsample, num_samples) and we must pad with sparse zeros
    log_probs = np.array(log_probs).T

    if downsampling:
        data = log_probs.flatten()
        row = np.repeat(downsample_indices, log_probs.shape[1])
        col = np.tile(list(range(log_probs.shape[1])), downsample)
        log_probs = coo_matrix((data, (row, col)), shape=(us.shape[0], len(unique_samples))) # (num_cells, num_samples)
        
    return log_probs
