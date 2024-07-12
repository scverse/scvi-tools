from collections.abc import Sequence

import torch
from anndata import AnnData
from torch.distributions import Distribution, Normal


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
):
    adata = self._validate_anndata(
        adata
    )  # need to maybe create class to inherit from vaemixin for now?
    dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

    qu_locs = []
    qu_scales = []
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

    qu_loc = torch.cat(qu_locs, 0).T
    qu_scale = torch.cat(
        qu_scales, 0
    ).T  # transpose because we need num cells to be rightmost dimension for mixture

    return Distribution.MixtureSameFamily(
        Distribution.Categorical(torch.ones(qu_loc.shape[1])), Normal(qu_loc, qu_scale)
    )


# def differential_abundance():
