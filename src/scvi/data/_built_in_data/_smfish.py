import logging
import os

import anndata
import pandas as pd

from scvi.utils import dependencies

logger = logging.getLogger(__name__)

_subtype_to_high_level_mapping = {
    "Astrocytes": ("Astrocyte Gfap", "Astrocyte Mfge8"),
    "Endothelials": ("Endothelial", "Endothelial 1"),
    "Inhibitory": (
        "Inhibitory Cnr1",
        "Inhibitory Kcnip2",
        "Inhibitory Pthlh",
        "Inhibitory Crhbp",
        "Inhibitory CP",
        "Inhibitory IC",
        "Inhibitory Vip",
    ),
    "Microglias": ("Perivascular Macrophages", "Microglia"),
    "Oligodendrocytes": (
        "Oligodendrocyte Precursor cells",
        "Oligodendrocyte COP",
        "Oligodendrocyte NF",
        "Oligodendrocyte Mature",
        "Oligodendrocyte MF",
    ),
    "Pyramidals": (
        "Pyramidal L2-3",
        "Pyramidal Cpne5",
        "Pyramidal L2-3 L5",
        "pyramidal L4",
        "Pyramidal L3-4",
        "Pyramidal Kcnip2",
        "Pyramidal L6",
        "Pyramidal L5",
        "Hippocampus",
    ),
}


@dependencies("pooch")
def _load_smfish(save_path: str = "data/", use_high_level_cluster=True) -> anndata.AnnData:
    import pooch

    save_path = os.path.abspath(save_path)
    adata = anndata.read_h5ad(
        pooch.retrieve(
            url="https://figshare.com/ndownloader/files/51096518",
            known_hash="a6bba682cf6804e4c1db07cbd2cb16a08143e0b814fd1bd1f936596aa1e27fd1",
            fname="smfish.h5ad",
            path=save_path,
            progressbar=True,
        )
    )
    if use_high_level_cluster:
        dataset = adata.obs.copy()
        dataset.str_labels = dataset.str_labels.astype(str)
        for high_level_cluster, subtypes in _subtype_to_high_level_mapping.items():
            dataset.loc[dataset.str_labels.isin(subtypes), "str_labels"] = high_level_cluster
        types_to_keep = [
            "Astrocytes",
            "Endothelials",
            "Inhibitory",
            "Microglias",
            "Oligodendrocytes",
            "Pyramidals",
        ]
        new_X = adata.X[dataset.str_labels.isin(types_to_keep)]
        dataset = dataset[dataset.str_labels.isin(types_to_keep)]
        dataset["str_labels"] = pd.Categorical(dataset["str_labels"])
        dataset["labels"] = pd.Categorical(dataset["str_labels"]).codes
        adata = anndata.AnnData(
            X=new_X,
            obs=dataset,
        )
        adata.var_names = adata.var_names
    adata.uns = {"cell_types": adata.obs.str_labels.cat.categories}
    return adata
