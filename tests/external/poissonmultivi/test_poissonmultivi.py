import numpy as np

from scvi.data import synthetic_iid
from scvi.external import POISSONMULTIVI


def test_poissonmultivi():
    data = synthetic_iid()
    POISSONMULTIVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = POISSONMULTIVI(
        data,
        n_genes=40,
        n_regions=60,
    )
    vae.train(1)
    vae.train(1, adversarial_mixing=False)
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_normalized_expression()
    vae.get_library_size_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.differential_expression(groupby="labels", group1="label_1")

    # Test with size factor expression
    data = synthetic_iid()
    data.obs["size_factor_expr"] = np.random.randint(1, 5, size=(data.shape[0],))
    POISSONMULTIVI.setup_anndata(data, batch_key="batch", size_factor_key_expr="size_factor_expr")
    vae = POISSONMULTIVI(
        data,
        n_genes=40,
        n_regions=60,
    )
    vae.train(3)

    # Test with size factor accessibility
    data = synthetic_iid()
    data.obs["size_factor_acc"] = np.random.randint(1, 5, size=(data.shape[0],))
    POISSONMULTIVI.setup_anndata(data, batch_key="batch", size_factor_key_acc="size_factor_acc")
    vae = POISSONMULTIVI(
        data,
        n_genes=40,
        n_regions=60,
    )
    vae.train(3)

    # Test with modality weights and penalties
    data = synthetic_iid()
    POISSONMULTIVI.setup_anndata(data, batch_key="batch")
    vae = POISSONMULTIVI(data, n_genes=40, n_regions=60, modality_weights="cell")
    vae.train(3)
    vae = POISSONMULTIVI(data, n_genes=40, n_regions=60, modality_weights="universal")
    vae.train(3)
    vae = POISSONMULTIVI(data, n_genes=40, n_regions=60, modality_penalty="MMD")
    vae.train(3)

    # Test with non-zero protein data
    data = synthetic_iid()
    POISSONMULTIVI.setup_anndata(
        data,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    vae = POISSONMULTIVI(
        data,
        n_genes=40,
        n_regions=60,
        modality_weights="cell",
    )
    assert vae.n_proteins == data.obsm["protein_expression"].shape[1]
    vae.train(3)
