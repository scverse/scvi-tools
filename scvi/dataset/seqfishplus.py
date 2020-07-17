import os
import zipfile
import anndata
import pandas as pd

from scvi.dataset import setup_anndata
from scvi.dataset._utils import _download


def seqfishplus(
    save_path="data/", tissue_region="subventricular cortex", run_setup_anndata=True
):

    if tissue_region == "subventricular cortex":
        file_prefix = "cortex_svz"
    elif tissue_region == "olfactory bulb":
        file_prefix = "ob"
    else:
        raise ValueError(
            '`tissue_type` must be "subventricular cortex" or "olfactory bulb", but got {}'.format(
                tissue_region
            )
        )

    save_path = os.path.abspath(save_path)
    url = "https://github.com/CaiGroup/seqFISH-PLUS/raw/master/sourcedata.zip"
    save_fn = "seqfishplus.zip"

    _download(url, save_path, save_fn)
    adata = _load_seqfishplus(
        os.path.join(save_path, save_fn), file_prefix, save_path, gene_by_cell=False
    )

    if run_setup_anndata:
        setup_anndata(adata)
    return adata


def _load_seqfishplus(path_to_file, file_prefix, save_path, gene_by_cell=False):
    counts_filename = "sourcedata/{}_counts.csv".format(file_prefix)
    coordinates_filename = "sourcedata/{}_cellcentroids.csv".format(file_prefix)
    extract_location = os.path.join(save_path, "seqfishplus")
    if not os.path.exists(extract_location):
        os.makedirs(extract_location)
    with zipfile.ZipFile(path_to_file) as f:
        f.extract(counts_filename, path=extract_location)
        f.extract(coordinates_filename, path=extract_location)

    df_counts = pd.read_csv(os.path.join(extract_location, counts_filename))
    adata = anndata.AnnData(df_counts)
    adata.var_names = df_counts.columns
    df_coordinates = pd.read_csv(os.path.join(extract_location, coordinates_filename))

    adata.obs["X"] = df_coordinates["X"].values
    adata.obs["Y"] = df_coordinates["Y"].values
    adata.obs["cell_id"] = df_coordinates["Cell ID"].values
    adata.obs["field_of_view"] = df_coordinates["Field of View"].values

    return adata
