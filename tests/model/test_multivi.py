import numpy as np

from scvi.data import synthetic_iid
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP
from scvi.model import (
    MULTIVI,
)

LEGACY_REGISTRY_KEYS = set(LEGACY_REGISTRY_KEY_MAP.values())
LEGACY_SETUP_DICT = {
    "scvi_version": "0.0.0",
    "categorical_mappings": {
        "_scvi_batch": {
            "original_key": "testbatch",
            "mapping": np.array(["batch_0", "batch_1"], dtype=object),
        },
        "_scvi_labels": {
            "original_key": "testlabels",
            "mapping": np.array(["label_0", "label_1", "label_2"], dtype=object),
        },
    },
    "extra_categoricals": {
        "mappings": {
            "cat1": np.array([0, 1, 2, 3, 4]),
            "cat2": np.array([0, 1, 2, 3, 4]),
        },
        "keys": ["cat1", "cat2"],
        "n_cats_per_key": [5, 5],
    },
    "extra_continuous_keys": np.array(["cont1", "cont2"], dtype=object),
    "data_registry": {
        "X": {"attr_name": "X", "attr_key": None},
        "batch_indices": {"attr_name": "obs", "attr_key": "_scvi_batch"},
        "labels": {"attr_name": "obs", "attr_key": "_scvi_labels"},
        "cat_covs": {
            "attr_name": "obsm",
            "attr_key": "_scvi_extra_categoricals",
        },
        "cont_covs": {
            "attr_name": "obsm",
            "attr_key": "_scvi_extra_continuous",
        },
    },
    "summary_stats": {
        "n_batch": 2,
        "n_cells": 400,
        "n_vars": 100,
        "n_labels": 3,
        "n_proteins": 0,
        "n_continuous_covs": 2,
    },
}


def test_multivi():
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(1, save_best=False)
    vae.train(1, adversarial_mixing=False)
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_normalized_expression()
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.differential_expression(groupby="labels", group1="label_1")

    # Test with size factor
    data = synthetic_iid()
    data.obs["size_factor"] = np.random.randint(1, 5, size=(data.shape[0],))
    MULTIVI.setup_anndata(data, batch_key="batch", size_factor_key="size_factor")
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(3)

    # Test with modality weights and penalties
    data = synthetic_iid()
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="cell")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="universal")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_penalty="MMD")
    vae.train(3)

    # Test with non-zero protein data
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
        modality_weights="cell",
    )
    assert vae.n_proteins == data.obsm["protein_expression"].shape[1]
    vae.train(3)
