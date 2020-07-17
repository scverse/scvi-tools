import logging
import pdb
import anndata
import os
import pandas as pd

from scvi.dataset._utils import _download
from scvi.dataset import setup_anndata

logger = logging.getLogger(__name__)


def seqfish(save_path="data/"):
    save_path = os.path.abspath(save_path)
    url = "https://www.cell.com/cms/attachment/2080562255/2072099886/mmc6.xlsx"
    save_fn = "SeqFISH.xlsx"
    _download(url, save_path, save_fn)
    adata = _load_seqfish_data(os.path.join(save_path, save_fn))
    return adata


def _load_seqfish_data(path_to_file):
    logger.info("Loading seqfish dataset from {}".format(path_to_file))
    xl = pd.ExcelFile(path_to_file)
    counts = xl.parse("Hippocampus Counts")
    X = counts.values[:, 1:].astype(int).T  # transpose because counts is genes X cells
    gene_names = counts.values[:, 0].astype(str)
    adata = anndata.AnnData(pd.DataFrame(data=X, columns=gene_names))
    logger.info("Finished loading seqfish dataset")
    return adata
