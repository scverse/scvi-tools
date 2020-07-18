import logging
import os
import scanpy as sc

from scvi.dataset._utils import _download
from scvi.dataset import setup_anndata

logger = logging.getLogger(__name__)


def breast_cancer_dataset(save_path="data/", run_setup_anndata=True):
    save_path = os.path.abspath(save_path)
    url = "http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/2016/07/Layer2_BC_count_matrix-1.tsv"
    save_fn = "Layer2_BC_count_matrix-1.tsv"
    _download(url, save_path, save_fn)
    adata = _load_csv(
        os.path.join(save_path, save_fn), delimiter="\t", gene_by_cell=False
    )
    if run_setup_anndata:
        setup_anndata(adata)
    return adata


def mouse_ob_dataset(save_path="data/", run_setup_anndata=True):
    save_path = os.path.abspath(save_path)
    url = "http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/2016/07/Rep11_MOB_count_matrix-1.tsv"
    save_fn = "Rep11_MOB_count_matrix-1.tsv"
    _download(url, save_path, save_fn)
    adata = _load_csv(
        os.path.join(save_path, save_fn), delimiter="\t", gene_by_cell=False
    )
    if run_setup_anndata:
        setup_anndata(adata)
    return adata


def _load_csv(path_to_file, gene_by_cell=False, delimiter=",", first_column_names=None):
    logger.info("Loading dataset from {}".format(path_to_file))
    adata = sc.read_csv(
        path_to_file, delimiter=delimiter, first_column_names=first_column_names
    )
    if gene_by_cell:
        adata.X = adata.X.T
    logger.info("Finished loading dataset")
    return adata
