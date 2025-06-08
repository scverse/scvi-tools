import logging
import os

from anndata import AnnData, read_h5ad

from scvi.utils import dependencies

logger = logging.getLogger(__name__)


@dependencies("pooch")
def _load_retina(save_path: str = "data/") -> AnnData:
    """Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells
    and 13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from
    the author. We also extract their normalized data with Combat and use it for benchmarking.
    """
    import pooch

    save_path = os.path.abspath(save_path)
    adata = read_h5ad(
        pooch.retrieve(
            url="https://figshare.com/ndownloader/files/51086201",
            known_hash="5363642ff02647d6868494b962ec962a5d2e3d90703415e245e7c1727c66cf21",
            fname="retina.h5ad",
            path=save_path,
            progressbar=True,
        )
    )
    return adata


@dependencies("pooch")
def _load_prefrontalcortex_starmap(save_path: str = "data/") -> AnnData:
    """Loads a starMAP dataset from the mouse pre-frontal cortex :cite:p:`Wang18`.

    Contains 3,704 cells and 166 genes.
    """
    import pooch

    save_path = os.path.abspath(save_path)
    adata = read_h5ad(
        pooch.retrieve(
            url="https://figshare.com/ndownloader/files/51086180",
            known_hash="c583eaef3835960405c6f1124f5fda36da80db3f940b76c9b2432a8d2e0b80ce",
            fname="mpfc-starmap.h5ad",
            path=save_path,
            progressbar=True,
        )
    )
    return adata


@dependencies("pooch")
def _load_frontalcortex_dropseq(save_path: str = "data/") -> AnnData:
    import pooch

    save_path = os.path.abspath(save_path)
    adata = read_h5ad(
        pooch.retrieve(
            url="https://figshare.com/ndownloader/files/51086207",
            known_hash="934a7179624a4c7c7dec1d5d53de5367fcd0054e5f19b7e245ecf2ecc88c188c",
            fname="fc-dropseq.h5ad",
            path=save_path,
            progressbar=True,
        )
    )
    # reorder labels such that layers of the cortex are in order
    # order_labels = [5, 6, 3, 2, 4, 0, 1, 8, 7, 9, 10, 11, 12, 13]
    # self.reorder_cell_types(self.cell_types[order_labels])
    return adata


@dependencies("pooch")
def _load_annotation_simulation(name: str, save_path: str = "data/") -> AnnData:
    """Simulated datasets for scANVI tutorials.

    Parameters
    ----------
    name
        One of the following:
        * ``'1'``
        * ``'2'``
        * ``'3'``
    save_path
        Location for saving the dataset.
    """
    import pooch

    if name == "1":
        fileid = "51086192"
        known_hash = "5d604adce93b3034885646605c2e9a72f5ccf8163caffb2930485f93a9fcb3a3"
    elif name == "2":
        fileid = "51086195"
        known_hash = "fdc2fb7c78e4c2a32877eb22aaed7cc627e22b256f122be670188c1069f741fa"
    else:
        fileid = "51086189"
        known_hash = "58c11e8c4134175c3f525f0d823a12420493cdf545f3904e0f09bec479c31e55"
    save_path = os.path.abspath(save_path)
    adata = read_h5ad(
        pooch.retrieve(
            url="https://figshare.com/ndownloader/files/" + fileid,
            known_hash=known_hash,
            fname=f"simulation_{name}.h5ad",
            path=save_path,
            progressbar=True,
        )
    )
    return adata
