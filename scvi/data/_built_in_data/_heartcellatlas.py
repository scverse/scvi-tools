import os

import anndata

from scvi.data import setup_anndata
from scvi.data._built_in_data._download import _download


def _load_heart_cell_atlas_subsampled(
    save_path: str = "data/",
    run_setup_anndata: bool = True,
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

    Returns
    -------
    AnnData
    """
    url = "https://github.com/YosefLab/scVI-data/raw/master/pbmc_10k_protein_v3.h5ad?raw=true"
    save_fn = "pbmc_10k_protein_v3.h5ad"
    _download(url, save_path, save_fn)
    dataset = anndata.read_h5ad(os.path.join(save_path, save_fn))

    if run_setup_anndata:
        setup_anndata(
            dataset,
            batch_key="batch",
        )

    return dataset
