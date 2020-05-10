import logging
import operator
import os
from functools import reduce

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from typing import Dict, Optional

from scvi.dataset.dataset import (
    DownloadableDataset,
    GeneExpressionDataset,
    CellMeasurement,
)

logger = logging.getLogger(__name__)


class AnnDatasetFromAnnData(GeneExpressionDataset):
    """Forms a ``GeneExpressionDataset`` from a ``anndata.AnnData`` object.

    :param ad: ``anndata.AnnData`` instance.
    :param batch_label: ``str`` representing AnnData obs column name for batches
    :param ctype_label: ``str`` representing AnnData obs column name for cell_types
    :param class_label: ``str`` representing AnnData obs column name for labels
    :param use_raw: if True, copies data from .raw attribute of AnnData
    """

    def __init__(
        self,
        ad: anndata.AnnData,
        batch_label: str = "batch_indices",
        ctype_label: str = "cell_types",
        class_label: str = "labels",
        use_raw: bool = False,
        cell_measurements_col_mappings: Optional[Dict[str, str]] = None,
    ):
        super().__init__()
        (
            X,
            batch_indices,
            labels,
            gene_names,
            cell_types,
            obs,
            obsm,
            var,
            _,
            uns,
        ) = extract_data_from_anndata(
            ad,
            batch_label=batch_label,
            ctype_label=ctype_label,
            class_label=class_label,
            use_raw=use_raw,
        )

        # Dataset API takes a dict as input
        obs = obs.to_dict(orient="list")
        var = var.to_dict(orient="list")

        # add external cell measurements
        Ys = []
        if cell_measurements_col_mappings is not None:
            for name, attr_name in cell_measurements_col_mappings.items():
                columns = uns[attr_name]
                measurement = CellMeasurement(
                    name=name,
                    data=obsm[name],
                    columns_attr_name=attr_name,
                    columns=columns,
                )
                Ys.append(measurement)

        self.populate_from_data(
            X=X,
            Ys=Ys,
            labels=labels,
            batch_indices=batch_indices,
            gene_names=gene_names,
            cell_types=cell_types,
            cell_attributes_dict=obs,
            gene_attributes_dict=var,
        )
        self.filter_cells_by_count()


class DownloadableAnnDataset(DownloadableDataset):
    """Forms a ``DownloadableDataset`` from a `.h5ad` file using the ``anndata`` package.

    :param filename: Name of the `.h5ad` file to save/load.
    :param save_path: Location to use when saving/loading the data.
    :param url: URL pointing to the data which will be downloaded
        if it's not already in ``save_path``.
    :param delayed_populating: Switch for delayed populating mechanism.
    :param batch_label: ``str`` representing AnnData obs column name for batches
    :param ctype_label: ``str`` representing AnnData obs column name for cell_types
    :param class_label: ``str`` representing AnnData obs column name for labels
    :param use_raw: if True, copies data from .raw attribute of AnnData

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
        batch_label: str = "batch_indices",
        ctype_label: str = "cell_types",
        class_label: str = "labels",
        use_raw: bool = False,
        cell_measurements_col_mappings: Optional[Dict[str, str]] = None,
    ):
        self.batch_label = batch_label
        self.ctype_label = ctype_label
        self.class_label = class_label
        self.use_raw = use_raw
        self.cell_measurements_col_mappings_temp = cell_measurements_col_mappings
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
            obs,
            obsm,
            var,
            _,
            uns,
        ) = extract_data_from_anndata(
            ad,
            batch_label=self.batch_label,
            ctype_label=self.ctype_label,
            class_label=self.class_label,
            use_raw=self.use_raw,
        )
        # Dataset API takes a dict as input
        obs = obs.to_dict(orient="list")
        var = var.to_dict(orient="list")

        # add external cell measurements
        Ys = []
        if self.cell_measurements_col_mappings_temp is not None:
            for name, attr_name in self.cell_measurements_col_mappings_temp.items():
                columns = uns[attr_name]
                measurement = CellMeasurement(
                    name=name,
                    data=obsm[name],
                    columns_attr_name=attr_name,
                    columns=columns,
                )
                Ys.append(measurement)

        self.populate_from_data(
            X=X,
            Ys=Ys,
            labels=labels,
            batch_indices=batch_indices,
            gene_names=gene_names,
            cell_types=cell_types,
            cell_attributes_dict=obs,
            gene_attributes_dict=var,
        )
        self.filter_cells_by_count()

        del self.cell_measurements_col_mappings_temp


def extract_data_from_anndata(
    ad: anndata.AnnData,
    batch_label: str = "batch_indices",
    ctype_label: str = "cell_types",
    class_label: str = "labels",
    use_raw: bool = False,
):
    data, labels, batch_indices, gene_names, cell_types = None, None, None, None, None
    # We use obs that will contain all the observation except those associated with
    #  batch_label, ctype_label and class_label.
    obs = ad.obs.copy()

    if use_raw:
        counts = ad.raw.X
    else:
        counts = ad.X

    # treat all possible cases according to anndata doc
    if isinstance(counts, np.ndarray):
        data = counts.copy()
    if isinstance(counts, pd.DataFrame):
        data = counts.values.copy()
    if sp_sparse.issparse(counts):
        # keep sparsity above 1 Gb in dense form
        if reduce(operator.mul, counts.shape) * counts.dtype.itemsize < 1e9:
            logger.info("Dense size under 1Gb, casting to dense format (np.ndarray).")
            data = counts.toarray()
        else:
            data = counts.copy()

    gene_names = np.asarray(ad.var.index.values, dtype=str)

    if batch_label in obs.columns:
        batch_indices = obs.pop(batch_label).values

    if ctype_label in obs.columns:
        cell_types = obs.pop(ctype_label)
        res = pd.factorize(cell_types)
        labels = res[0].astype(int)
        cell_types = np.array(res[1]).astype(str)

    elif class_label in obs.columns:
        labels = obs.pop(class_label)

    return (
        data,
        batch_indices,
        labels,
        gene_names,
        cell_types,
        obs,
        ad.obsm,
        ad.var,
        ad.varm,
        ad.uns,
    )
