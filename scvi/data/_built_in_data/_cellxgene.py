import os
from typing import Optional

import pandas as pd
import requests
from anndata import AnnData, read_h5ad

from scvi.data._download import _download

CELLXGENE_PRODUCTION_HOST = "api.cellxgene.cziscience.com"
CELLXGENE_PRODUCTION_ENDPOINT = f"https://{CELLXGENE_PRODUCTION_HOST}"
DATASETS = f"{CELLXGENE_PRODUCTION_ENDPOINT}/dp/v1/datasets/"
COLLECTIONS = f"{CELLXGENE_PRODUCTION_ENDPOINT}/dp/v1/collections/"
CELLXGENE_EXPLORER = "https://cellxgene.cziscience.com/e/"


def _presign_url(url: str):
    resp = requests.post(url)
    return resp.json()["presigned_url"]


def _load_cellxgene_dataset(
    url: str,
    filename: Optional[str] = None,
    save_path: str = "data/",
) -> AnnData:
    """
    Loads a file from `cellxgene <https://cellxgene.cziscience.com/>`_ portal.

    Parameters
    ----------
    url
        URL to cellxgene session
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
    dataset_id = url.split("/")[-1][:-4]
    file_id = None

    # This is a bit of an inefficient way to do this, but it's the only way to get the file id
    # from the cellxgene explorer link
    # We have to search through all the collections to find the file id
    collections_json = requests.get(COLLECTIONS).json()
    db_tbl = pd.DataFrame.from_records(collections_json["collections"])
    for _, rec in db_tbl.iterrows():
        rec_resp = requests.get(COLLECTIONS + rec["id"]).json()
        datasets = rec_resp["datasets"]
        collection_assets = []
        for i in datasets:
            for j in i["dataset_assets"]:
                if j["filetype"] == "H5AD":
                    collection_assets.append(j)
        # see if the dataset is in this collection
        for asset in collection_assets:
            if asset["dataset_id"] == dataset_id:
                file_id = asset["id"]
                break
    if file_id is None:
        raise ValueError("Dataset not found in any collection")

    url = f"{DATASETS}{dataset_id}/asset/{file_id}"
    presigned_url = _presign_url(url)
    _download(presigned_url, save_path, filename)
    file_path = os.path.join(save_path, filename)
    adata = read_h5ad(file_path)

    return adata
