import os
import re
from typing import Optional, Union

from anndata import AnnData, read_h5ad

from scvi._decorators import dependencies


def _parse_dataset_id(url: str):
    match = re.search(r"/e/(.+)", url)
    if match:
        return match.group(1).split(".")[0]
    else:
        raise ValueError(f"Could not parse dataset id from url {url}")


@dependencies("cellxgene_census")
def _load_cellxgene_dataset(
    url: str,
    filename: Optional[str] = None,
    save_path: str = "data/",
    return_path: bool = False,
) -> Union[AnnData, str]:
    """Loads a file from `cellxgene <https://cellxgene.cziscience.com/>`_ portal.

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
    return_path
        If True, returns the path to the downloaded file instead of the AnnData object.

    Returns
    -------
    adata initialized with cellxgene data
    """
    from cellxgene_census import download_source_h5ad

    dataset_id = _parse_dataset_id(url)

    if filename is None:
        filename = f"{dataset_id}.h5ad"
    file_path = os.path.join(save_path, filename)
    # check if file exists
    if not os.path.exists(file_path):
        download_source_h5ad(dataset_id, to_path=file_path)
    if return_path:
        return file_path
    adata = read_h5ad(file_path)
    return adata
