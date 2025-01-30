from collections.abc import Sequence

import numpy as np

import pandas as pd
import torch
import torch.distributions as dist
from anndata import AnnData
from torch import Tensor
from tqdm import tqdm


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    sample: str | int | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
    dof: float | None = 3.,
    qu_scale = None,
    qu_loc = None,
) -> dist.Distribution:
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

    qu_loc = torch.tensor(qu_loc, device='cuda').T
    qu_scale = torch.tensor(qu_scale, device='cuda').T
    
    if dof is None:
        components = dist.Normal(qu_loc, qu_scale)
    else:
        components = dist.StudentT(dof, qu_loc, qu_scale)
    return dist.MixtureSameFamily(
        dist.Categorical(logits=torch.ones(qu_loc.shape[1], device='cuda')), components)

def differential_abundance(
    self,
    adata: AnnData | None = None,
    sample_key: str | None = None,
    batch_size: int = 128,
    downsample_cells: int | None = None,
    dof: float | None = None,
    qu_loc = None,
    qu_scale = None,
) -> pd.DataFrame:
    adata = self._validate_anndata(adata)

    us = self.get_latent_representation(
        batch_size=batch_size, return_dist=False, give_mean=True
    )
    
    unique_samples = adata.obs[sample_key].unique()
    dataloader = torch.utils.data.DataLoader(us, batch_size=batch_size)
    log_probs = []
    for sample_name in tqdm(unique_samples):
        indices = np.where(adata.obs[sample_key] == sample_name)[0]
        if downsample_cells is not None and downsample_cells < indices.shape[0]:
            indices = np.random.choice(indices, downsample_cells, replace=False)

        ap = get_aggregated_posterior(self, adata=adata, indices=indices, dof=dof, qu_loc=qu_loc, qu_scale=qu_scale)
        log_probs_ = []
        for u_rep in dataloader:
            u_rep = u_rep.to('cuda')
            log_probs_.append(ap.log_prob(u_rep).sum(-1, keepdims=True))
        log_probs.append(torch.cat(log_probs_, axis=0).cpu().numpy())

    log_probs = np.concatenate(log_probs, 1)
    log_probs_df = pd.DataFrame(data=log_probs, index=adata.obs_names.to_numpy(), columns=unique_samples)
    return log_probs_df
