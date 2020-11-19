import logging
from typing import Optional

import anndata
import numpy as np
import torch

from scvi import _CONSTANTS
from scvi.core.modules import TOTALVAE

from .scvi_data_loader import ScviDataLoader

logger = logging.getLogger(__name__)


def _unpack_tensors(tensors):
    x = tensors[_CONSTANTS.X_KEY]
    local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
    local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
    batch_index = tensors[_CONSTANTS.BATCH_KEY]
    labels = tensors[_CONSTANTS.LABELS_KEY]
    y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
    return x, local_l_mean, local_l_var, batch_index, labels, y


class TotalDataLoader(ScviDataLoader):
    """
    Extended data loader for totalVI.

    Parameters
    ----------
    model :
        A model instance from class ``TOTALVI``
    adata:
        A registered AnnData object
    shuffle :
        Specifies if a `RandomSampler` or a `SequentialSampler` should be used
    indices :
        Specifies how the data should be split with regards to train/test or labelled/unlabelled
    use_cuda :
        Default: ``True``
    data_loader_kwargs :
        Keyword arguments to passed into the `DataLoader`

    """

    def __init__(
        self,
        model: TOTALVAE,
        adata: anndata.AnnData,
        shuffle: bool = False,
        indices: Optional[np.ndarray] = None,
        use_cuda: bool = True,
        batch_size: int = 256,
        data_loader_kwargs=dict(),
    ):
        super().__init__(
            model,
            adata,
            shuffle=shuffle,
            indices=indices,
            use_cuda=use_cuda,
            batch_size=batch_size,
            data_loader_kwargs=data_loader_kwargs,
        )

    @property
    def _data_and_attributes(self):
        return {
            _CONSTANTS.X_KEY: np.float32,
            _CONSTANTS.BATCH_KEY: np.int64,
            _CONSTANTS.LOCAL_L_MEAN_KEY: np.float32,
            _CONSTANTS.LOCAL_L_VAR_KEY: np.float32,
            _CONSTANTS.LABELS_KEY: np.int64,
            _CONSTANTS.PROTEIN_EXP_KEY: np.float32,
        }

    # @torch.no_grad()
    # def elbo(self):
    #     elbo = self.compute_elbo()
    #     return elbo

    # elbo.mode = "min"

    @torch.no_grad()
    def reconstruction_error(self, mode="total"):
        ll_gene, ll_protein = self.compute_reconstruction_error(self.model)
        if mode == "total":
            return ll_gene + ll_protein
        elif mode == "gene":
            return ll_gene
        else:
            return ll_protein

    reconstruction_error.mode = "min"

    @torch.no_grad()
    def marginal_ll(self, n_mc_samples=1000):
        ll = self.compute_marginal_log_likelihood()
        return ll

    @torch.no_grad()
    def get_protein_background_mean(self):
        background_mean = []
        for tensors in self:
            x, _, _, batch_index, label, y = _unpack_tensors(tensors)
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=label, n_samples=1
            )
            b_mean = outputs["py_"]["rate_back"]
            background_mean += [np.array(b_mean.cpu())]
        return np.concatenate(background_mean)
