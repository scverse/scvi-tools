import os
from typing import Optional, Union

import requests
from anndata import AnnData, read_h5ad

from scvi.data._download import _download

BACKEND = "https://public-backend.production.single-cell.czi.technology/curation/v1/collections/"


def _load_cellxgene_dataset(
    url: str,
    collection_id: Optional[str] = None,
    filename: Optional[str] = None,
    save_path: str = "data/",
    return_path: bool = False,
) -> Union[AnnData, str]:
    """
    Loads a file from `cellxgene <https://cellxgene.cziscience.com/>`_ portal.

    Parameters
    ----------
    url
        URL to cellxgene explorer
    collection_id
        Cellxgene collection ID. If None, will find the correct collection.
    filename
        manual override of the filename to write to.
    save_path
        Location to use when saving/loading the data.

    Returns
    -------
    adata initialized with cellxgene data

    Notes
    -----
    API here:

    https://public-backend.production.single-cell.czi.technology/
    curation/ui/#/collection/backend.portal.api.curation.v1.curation.collections.
    collection_id.datasets.dataset_id.assets.get
    """
    # get the dataset id from the url and remove .cxg
    split_url = url.split("/")
    dataset_id = split_url[-2] if split_url[-1] == "" else split_url[-1]
    dataset_id = dataset_id.split(".")[0]

    # the get request actually ignores the collection id as dataset ids are unique
    collection_id = "random" if collection_id is None else collection_id
    rec = requests.get(f"{BACKEND}/{collection_id}/datasets/{dataset_id}/assets")
    rec.raise_for_status()
    assets = rec.json()
    presigned_url = None
    for asset in assets:
        if asset["filename"].endswith(".h5ad"):
            presigned_url = asset["presigned_url"]
    if presigned_url is None:
        raise ValueError("No h5ad file found in dataset")
    if filename is None:
        filename = "local.h5ad"
    _download(presigned_url, save_path, filename)
    file_path = os.path.join(save_path, filename)
    if return_path:
        return file_path
    adata = read_h5ad(file_path)
    return adata
