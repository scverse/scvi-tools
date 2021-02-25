import os

import anndata

from scvi.data import setup_anndata
from scvi.data._built_in_data._download import _download


def _load_heart_cell_atlas_subsampled(
    save_path: str = "data/",
    run_setup_anndata: bool = True,
    remove_nuisance_clusters: bool = True,
):
    """
    Combined single cell and single nuclei RNA-Seq data of 485K cardiac cells with annotations.

    Dataset was filtered down randomly to 20k cells using :func:`~scanpy.pp.subsample`. The original
    data can be sourced from https://www.heartcellatlas.org/#DataSources.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning
    remove_nuisance_clusters
        Remove doublets and unsassigned cells

    Returns
    -------
    AnnData

    Notes
    -----
    The data were filtered using the following sequence::

        >>> adata = anndata.read_h5ad(path_to_anndata)
        >>> bdata = sc.pp.subsample(adata, n_obs=20000, copy=True)
        >>> sc.pp.filter_genes(bdata, min_counts=3)
        >>> bdata.write_h5ad(path, compression="gzip")
    """
    url = "https://github.com/YosefLab/scVI-data/blob/master/hca_subsampled_20k.h5ad?raw=true"
    save_fn = "hca_subsampled_20k.h5ad"
    _download(url, save_path, save_fn)
    dataset = anndata.read_h5ad(os.path.join(save_path, save_fn))

    if remove_nuisance_clusters:
        remove = ["doublets", "NotAssigned"]
        keep = [c not in remove for c in dataset.obs.cell_type.values]
        dataset = dataset[keep, :].copy()

    if run_setup_anndata:
        setup_anndata(
            dataset,
        )

    return dataset
