import logging
import operator
import os

from functools import reduce
from typing import List

import anndata
import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix

from scvi.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class AnnDataset(DownloadableDataset):
    r"""Loads a `.h5ad` file using the ``anndata`` package.
    Also supports loading an ``AnnData`` object directly.

    Args:
        :filename_or_anndata: Name of the `.h5ad` file or an existing ``AnnData``.
        :save_path: Save path of the dataset. Default: ``'data/'``.
        :url: Url of the remote dataset. Default: ``None``.
        :new_n_genes: Number of subsampled genes. Default: ``False``.
        :subset_genes: List of genes for subsampling. Default: ``None``.


    Examples:
        >>> # Loading a local dataset
        >>> local_ann_dataset = AnnDataset("TM_droplet_mat.h5ad", save_path = 'data/')

    .. _Anndata:
        http://anndata.readthedocs.io/en/latest/

    """

    def __init__(
        self,
        filename: str,
        save_path: str = "data/",
        url: str = None,
        delayed_populating: bool = False,
        preprocess: bool = True,
    ):
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
            preprocess=preprocess,
        )

    def load_from_disk(self):
        ad = anndata.read_h5ad(
            os.path.join(self.save_path, self.filenames[0])
        )  # obs = cells, var = genes

        # extract GeneExpressionDataset relevant attributes
        # and provide access to annotations from the underlying AnnData object.
        (X,
         batch_indices,
         labels,
         gene_names,
         cell_types,
         self.obs,
         self.obsm,
         self.var,
         self.varm
         ) = extract_data_from_anndata(ad)
        return {
            "X": X,
            "batch_indices": batch_indices,
            "labels": labels,
            "gene_names": gene_names,
            "cell_types": cell_types,
        }

    def preprocess(self, **kwargs):
        return kwargs

    def instantiate_gene_expression_dataset(self, **kwargs):
        super(DownloadableDataset, self).__init__(**kwargs)

    @classmethod
    def from_anndata(cls, ad: anndata.AnnData,):
        data, gene_names, batch_indices, cell_types, labels, obs, obsm, var, varm = extract_data_from_anndata(ad)

        dataset = super(DownloadableDataset, cls).__init__(
            X=data,
            batch_indices=batch_indices,
            labels=labels,
            cell_types=cell_types,
            gene_names=gene_names,

        )

        dataset.obs = obs
        dataset.obsm = obsm
        dataset.var = var
        dataset.varm = varm

        return dataset


def extract_data_from_anndata(ad: anndata.AnnData):
    data, gene_names, batch_indices, cell_types, labels = None, None, None, None, None

    # treat all possible cases according to anndata doc
    if isinstance(ad.X, np.ndarray):
        data = ad.X.copy()
    if isinstance(ad.X, pd.DataFrame):
        data = ad.X.values
    if isinstance(ad.X, csr_matrix):
        # keep sparsity above 1 Gb in dense form
        if reduce(operator.mul, ad.X.shape) * ad.X.dtype.itemsize < 1e9:
            data = ad.X.toarray()
        else:
            data = ad.X.copy()

    gene_names = np.array(ad.var.index.values, dtype=str)

    if 'batch_indices' in ad.obs.columns:
        batch_indices = ad.obs['batch_indices'].values

    if 'cell_types' in ad.obs.columns:
        cell_types = ad.obs['cell_types']
        cell_types = cell_types.drop_duplicates().values.astype('str')

    if 'labels' in ad.obs.columns:
        labels = ad.obs['labels']
    elif 'cell_types' in ad.obs.columns:
        labels = ad.obs['cell_types'].rank(method='dense').values.astype('int')

    return data, batch_indices, labels, gene_names, cell_types, ad.obs, ad.obsm, ad.var, ad.varm
