import os
from pathlib import Path
from typing import Union

import pandas as pd
from anndata import AnnData
from scipy.io import mmread


def read_10x_atac(base_path: Union[str, Path]) -> AnnData:
    """
    Read scATAC-seq data outputted by 10x Genomics software.

    Parameters
    ----------
    base_path
        Path to directory with matrix, bed file, etc.
    """
    data = mmread(os.path.join(base_path, "matrix.mtx")).transpose()
    coords = pd.read_csv(
        os.path.join(base_path, "peaks.bed"),
        sep="\t",
        header=None,
        index_col=None,
    )
    coords.rename({0: "chr", 1: "start", 2: "end"}, axis="columns", inplace=True)
    coords.set_index(
        coords.chr.astype(str)
        + ":"
        + coords.start.astype(str)
        + "-"
        + coords.end.astype(str),
        inplace=True,
    )
    coords.index = coords.index.astype(str)

    cell_annot = pd.read_csv(
        os.path.join(base_path, "barcodes.tsv"), sep="-", header=None, index_col=None
    )
    cell_annot.rename({0: "barcode", 1: "batch_id"}, axis="columns", inplace=True)
    cell_annot.set_index("barcode", inplace=True)
    cell_annot.index = cell_annot.index.astype(str)

    return AnnData(data.tocsr(), var=coords, obs=cell_annot)


def read_10x_multiome(base_path: Union[str, Path]) -> AnnData:
    """
    Read Multiome (scRNA + scATAC) data outputted by 10x Genomics software.

    Parameters
    ----------
    base_path
        Path to directory with matrix, barcodes file, etc.
    """
    data = mmread(os.path.join(base_path, "matrix.mtx")).transpose()

    features = pd.read_csv(
        os.path.join(base_path, "features.tsv"),
        sep="\t",
        header=None,
        index_col=1,
    )
    features.rename(
        {0: "ID", 2: "modality", 3: "chr", 4: "start", 5: "end"},
        axis="columns",
        inplace=True,
    )
    features.index.name = None

    cell_annot = pd.read_csv(
        os.path.join(base_path, "barcodes.tsv"), sep="-", header=None, index_col=None
    )
    cell_annot.rename({0: "barcode", 1: "batch_id"}, axis="columns", inplace=True)
    cell_annot.set_index("barcode", inplace=True)
    cell_annot.index = cell_annot.index.astype(str)

    return AnnData(data.tocsr(), var=features, obs=cell_annot)
