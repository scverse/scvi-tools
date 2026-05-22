import itertools

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from mudata import MuData
from scipy import sparse as sp_sparse

from scvi.data import synthetic_iid


# helper function for testing purposes ; could be moved in the same file as generate_synthetic()
def _generate_synthetic_hierarchy(
    batch_size: int = 128,
    n_genes: int = 100,
    n_proteins: int = 100,
    n_batches: int = 2,
    n_labels_1: int = 2,
    n_labels_2: int = 11,
    n_sites: int = 2,
    sparse: bool = False,
) -> AnnData:
    """
    New method to generate test data with two-layer labels.
    """
    n_total = batch_size * n_batches * n_sites
    data = np.random.negative_binomial(5, 0.3, size=(n_total))
    mask = np.random.binomial(n=1, p=0.7, size=(n_total, n_genes))
    data = data * mask  # We put the batch index first
    labels_1 = np.random.randint(0, n_labels_1, size=(n_total,))
    labels_1 = np.array([f"label_{i}" for i in labels_1])
    labels_2 = np.random.randint(0, n_labels_2 // n_sites, size=(n_total,))

    batch = []
    site = []
    labels_2 = []
    for site, batch in itertools.product(range(n_sites), range(n_batches)):
        batch += [f"batch_{batch}_site_{site}"] * batch_size
        site += [f"site_{site}"] * batch_size
    labels_1 = np.array([f"label_{d}_{s}" for d, s in zip(labels_2, site, strict=True)])

    if sparse:
        data = sp_sparse.csr_matrix(data)
    adata = AnnData(data)
    adata.obs["batch"] = pd.Categorical(batch)
    adata.obs["site"] = pd.Categorical(site)
    adata.obs["labels_1"] = pd.Categorical(labels_1)
    adata.obs["labels_2"] = pd.Categorical(labels_2)

    # Protein measurements
    p_data = np.random.negative_binomial(5, 0.3, size=(adata.shape[0], n_proteins))
    adata.obsm["protein_expression"] = p_data
    adata.uns["protein_names"] = np.arange(n_proteins).astype(str)

    return adata


# helper function for testing purposes ; could be moved int he same file as synthetic_iid()
def synthetic_iid_hierarchy(
    batch_size: int | None = 200,
    n_genes: int | None = 100,
    n_proteins: int | None = 100,
    n_batches: int | None = 2,
    n_sites: int | None = 2,
    n_labels_1: int | None = 2,
    n_labels_2: int | None = 10,
    sparse: bool = False,
) -> AnnData:
    """Synthetic dataset with ZINB distributed RNA and NB distributed protein, with three-layer annotation.
    This dataset is just for testing purposed and not meant for modeling or research.
    Each value is independently and identically distributed.
    Parameters
    ----------
    batch_size
        Number of cells per batch
    n_genes
        Number of genes
    n_proteins
        Number of proteins
    n_batches
        Number of batches
    n_sites
        Number of sites
    n_labels_1
        Number of cell types of layer 1
    n_labels_2
        Number of cell types of layer 2
    sparse
        Whether to use a sparse matrix
    Returns
    -------
    AnnData with batch info (``.obs['batch']``), label info (``.obs['labels']``),
    site info (``.obs['site']``),
    protein expression (``.obsm["protein_expression"]``) and
    protein names (``.obs['protein_names']``)
    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.synthetic_iid()
    """

    return _generate_synthetic_hierarchy(
        batch_size=batch_size,
        n_genes=n_genes,
        n_proteins=n_proteins,
        n_batches=n_batches,
        n_sites=n_sites,
        n_labels_1=n_labels_1,
        n_labels_2=n_labels_2,
        sparse=sparse,
    )


def test_methylvi():
    adata1 = synthetic_iid()
    adata1.layers["mc"] = adata1.X
    adata1.layers["cov"] = adata1.layers["mc"] + 10

    adata2 = synthetic_iid()
    adata2.layers["mc"] = adata2.X
    adata2.layers["cov"] = adata2.layers["mc"] + 10

    mdata = MuData({"mod1": adata1, "mod2": adata2})

    METHYLVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        methylation_contexts=["mod1", "mod2"],
        batch_key="batch",
        modalities={"batch_key": "mod1"},
    )
    vae = METHYLVI(
        mdata,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_normalized_methylation()  # Retrieve methylation for all contexts
    vae.get_normalized_methylation(context="mod1")  # Retrieve for specific context
    with pytest.raises(ValueError):  # Should fail when invalid context selected
        vae.get_normalized_methylation(context="mod3")
    vae.get_latent_representation()
    vae.differential_methylation(groupby="mod1:labels", group1="label_1")
