import os
from pathlib import Path
from typing import Union

import pandas as pd
from anndata import AnnData
from scipy.io import mmread


def read_10x_atac(base_path: Union[str, Path]) -> AnnData:
    """
    Read scATAC-seq data ouputted by 10x Genomics software.

    Parameters
    ----------
    base_path
        Path to directory with matrix, bed file, etc.
    """
    data = mmread(os.path.join(base_path, "matrix.mtx")).transpose()
    coords = pd.read_csv(
        os.path.join(base_path, "peaks.bed"), sep="\t", header=None, index_col=None
    )
    coords.rename({0: "Chromosome", 1: "Start", 2: "End"}, axis="columns", inplace=True)
    cell_annot = pd.DataFrame(
        [
            line.strip().split("-")
            for line in open(os.path.join(base_path, "barcodes.tsv"), "r").readlines()
        ]
    )
    cell_annot.rename({0: "cell_barcode", 1: "batch_id"}, axis="columns", inplace=True)
    cell_annot = cell_annot.set_index("cell_barcode")
    return AnnData(data.tocsr(), var=coords, obs=cell_annot)
