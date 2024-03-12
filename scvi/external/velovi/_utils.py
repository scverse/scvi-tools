from __future__ import annotations

from pathlib import Path

import pandas as pd
from anndata import AnnData

from scvi.utils import dependencies


def get_permutation_scores(save_path: str | Path = Path("data/")) -> pd.DataFrame:
    """``BETA`` Get the reference permutation scores on positive and negative controls.

    Parameters
    ----------
    save_path
        Path to save the csv file.
    """
    from urllib.request import urlretrieve

    if isinstance(save_path, str):
        save_path = Path(save_path)
    save_path.mkdir(parents=True, exist_ok=True)

    if not (save_path / "permutation_scores.csv").is_file():
        URL = "https://figshare.com/ndownloader/files/36658185"
        urlretrieve(url=URL, filename=save_path / "permutation_scores.csv")

    return pd.read_csv(save_path / "permutation_scores.csv")


@dependencies("scvelo")
def preprocess_data(
    adata: AnnData,
    spliced_layer: str | None = "Ms",
    unspliced_layer: str | None = "Mu",
    min_max_scale: bool = True,
    filter_on_r2: bool = True,
) -> AnnData:
    """``BETA`` Preprocess data.

    This function removes poorly detected genes and minmax scales the data.

    Parameters
    ----------
    adata
        Annotated data matrix.
    spliced_layer
        Name of the spliced layer.
    unspliced_layer
        Name of the unspliced layer.
    min_max_scale
        Min-max scale spliced and unspliced.
    filter_on_r2
        Filter out genes according to linear regression fit.

    Returns
    -------
    Preprocessed adata.
    """
    import numpy as np
    import scvelo as scv
    from sklearn.preprocessing import MinMaxScaler

    if min_max_scale:
        scaler = MinMaxScaler()
        adata.layers[spliced_layer] = scaler.fit_transform(adata.layers[spliced_layer])

        scaler = MinMaxScaler()
        adata.layers[unspliced_layer] = scaler.fit_transform(adata.layers[unspliced_layer])

    if filter_on_r2:
        scv.tl.velocity(adata, mode="deterministic")

        adata = adata[
            :, np.logical_and(adata.var.velocity_r2 > 0, adata.var.velocity_gamma > 0)
        ].copy()
        adata = adata[:, adata.var.velocity_genes].copy()

    return adata
