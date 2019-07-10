import logging
import operator
import os
from functools import reduce

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import DownloadableDataset, GeneExpressionDataset

logger = logging.getLogger(__name__)


class AnnDatasetFromAnnData(GeneExpressionDataset):
    """Forms a ``GeneExpressionDataset`` from a ``anndata.AnnData`` object.

    :param ad: ``anndata.AnnData`` instance.
    """

    def __init__(self, ad: anndata.AnnData):
        super().__init__()
        (
            X,
            batch_indices,
            labels,
            gene_names,
            cell_types,
            self.obs,
            self.obsm,
            self.var,
            self.varm,
            self.uns
        ) = extract_data_from_anndata(ad)
        self.populate_from_data(
            X=X,
            batch_indices=batch_indices,
            gene_names=gene_names,
            cell_types=cell_types,
        )
        self.filter_cells_by_count()


class DownloadableAnnDataset(DownloadableDataset):
    """Forms a ``DownloadableDataset`` from a `.h5ad` file using the ``anndata`` package.

    :param filename: Name of the `.h5ad` file to save/load.
    :param save_path: Location to use when saving/loading the data.
    :param url: URL pointing to the data which will be downloaded
        if it's not already in ``save_path``.
    :param delayed_populating: Switch for delayed populating mechanism.

        Examples:
        >>> # Loading a local dataset
        >>> dataset = DownloadableAnnDataset("TM_droplet_mat.h5ad", save_path = 'data/')

    .. _Anndata:
        http://anndata.readthedocs.io/en/latest/
    """

    def __init__(
        self,
        filename: str = "anndataset",
        save_path: str = "data/",
        url: str = None,
        delayed_populating: bool = False,
    ):
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        ad = anndata.read_h5ad(
            os.path.join(self.save_path, self.filenames[0])
        )  # obs = cells, var = genes

        # extract GeneExpressionDataset relevant attributes
        # and provide access to annotations from the underlying AnnData object.
        (
            X,
            batch_indices,
            labels,
            gene_names,
            cell_types,
            self.obs,
            self.obsm,
            self.var,
            self.varm,
            self.uns,
        ) = extract_data_from_anndata(ad)
        self.populate_from_data(
            X=X,
            batch_indices=batch_indices,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
        )
        self.filter_cells_by_count()


def extract_data_from_anndata(ad: anndata.AnnData):
    data, labels, batch_indices, gene_names, cell_types = None, None, None, None, None

    # treat all possible cases according to anndata doc
    if isinstance(ad.X, np.ndarray):
        data = ad.X.copy()
    if isinstance(ad.X, pd.DataFrame):
        data = ad.X.values
    if isinstance(ad.X, csr_matrix):
        # keep sparsity above 1 Gb in dense form
        if reduce(operator.mul, ad.X.shape) * ad.X.dtype.itemsize < 1e9:
            logger.info("Dense size under 1Gb, casting to dense format (np.ndarray).")
            data = ad.X.toarray()
        else:
            data = ad.X.copy()

    gene_names = np.asarray(ad.var.index.values, dtype=str)

    if "batch_indices" in ad.obs.columns:
        batch_indices = ad.obs["batch_indices"].values

    if "cell_types" in ad.obs.columns:
        cell_types = ad.obs["cell_types"]
        labels = cell_types.rank(method="dense").values.astype("int")
        cell_types = cell_types.drop_duplicates().values.astype("str")

    if "labels" in ad.obs.columns:
        labels = ad.obs["labels"]

    return (
        data,
        batch_indices,
        labels,
        gene_names,
        cell_types,
        ad.obs,
        ad.obsm,
        ad.var,
        ad.varm,
        ad.uns,
    )
