import os

import numpy as np
import pytest
import scanpy as sc

from scvi.data import synthetic_iid
from scvi.model import CondSCVI, DestVI


def test_destvi():
    # Step1 learn CondSCVI
    n_latent = 2
    n_labels = 5
    n_layers = 2
    dataset = synthetic_iid(n_labels=n_labels)
    dataset.obs["overclustering_vamp"] = list(range(dataset.n_obs))
    CondSCVI.setup_anndata(dataset, labels_key="labels", batch_key="batch")
    sc_model = CondSCVI(
        dataset, n_latent=n_latent, n_layers=n_layers, prior="mog", num_classes_mog=10
    )
    sc_model.train(1, train_size=1)

    sc_model.get_normalized_expression(dataset)
    sc_model.get_elbo()
    sc_model.get_reconstruction_error()
    sc_model.get_latent_representation()
    sc_model.get_vamp_prior(dataset, p=100)

    # step 2 Check model setup
    DestVI.setup_anndata(dataset, layer=None)

    # Test clustering outside of get_vamp_prior

    _ = DestVI.from_rna_model(dataset, sc_model, vamp_prior_p=dataset.n_obs)

    del dataset.obs["overclustering_vamp"]

    # step 3 learn destVI with multiple amortization scheme

    for amor_scheme in ["both", "latent"]:
        DestVI.setup_anndata(dataset, layer=None)
        spatial_model = DestVI.from_rna_model(dataset, sc_model, amortization=amor_scheme)
        spatial_model.view_anndata_setup()
        spatial_model.train(max_epochs=1)
        assert not np.isnan(spatial_model.history["elbo_train"].values[0][0])

        assert spatial_model.get_proportions().shape == (dataset.n_obs, n_labels)
        assert spatial_model.get_gamma(return_numpy=True).shape == (
            dataset.n_obs,
            n_latent,
            n_labels,
        )

        assert spatial_model.get_scale_for_ct("label_0", np.arange(50)).shape == (
            50,
            dataset.n_vars,
        )

        with pytest.raises(NotImplementedError):
            spatial_model.get_normalized_expression()


@pytest.mark.internet
def test_destvi_new(save_path: str):
    CELL_TYPE_ID = "broad_cell_types"
    G = 2000
    CELL_TYPE_HIGHRES_ID = "cell_types"
    BATCH_KEY = "batch"
    MODELS_PATH = os.path.join(save_path, "models")
    VISIUM_COLON_SC_MODEL_PATH = os.path.join(MODELS_PATH, "Visium_DestVI_v2_sc")
    VISIUM_COLON_ST_MODEL_PATH = os.path.join(MODELS_PATH, "Visium_DestVI_v2_st")
    st_adata_path = os.path.join(save_path, "st_lymph_node_preprocessed.h5ad")
    st_adata = sc.read(
        st_adata_path,
        backup_url="https://exampledata.scverse.org/scvi-tools/st_lymph_node_preprocessed.h5ad",
    )

    sc_adata_path = os.path.join(save_path, "sc_lymph_node_preprocessed.h5ad")
    sc_adata = sc.read(
        sc_adata_path,
        backup_url="https://exampledata.scverse.org/scvi-tools/sc_lymph_node_preprocessed.h5ad",
    )

    # let us filter some genes
    sc.pp.filter_genes(sc_adata, min_counts=10)
    sc_adata.layers["counts"] = sc_adata.X.copy()
    sc.pp.highly_variable_genes(
        sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3"
    )
    sc.pp.normalize_total(sc_adata, target_sum=10e4)
    sc.pp.log1p(sc_adata)
    sc_adata.raw = sc_adata

    st_adata.layers["counts"] = st_adata.X.copy()
    st_adata.obsm["spatial"] = st_adata.obsm["location"]
    sc.pp.normalize_total(st_adata, target_sum=10e4)
    sc.pp.log1p(st_adata)
    st_adata.raw = st_adata

    # filter genes to be the same on the spatial data
    intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
    st_adata = st_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()
    G = len(intersect)
    print(G)

    CondSCVI.setup_anndata(
        sc_adata,
        layer="counts",
        labels_key=CELL_TYPE_ID,
        fine_labels_key=CELL_TYPE_HIGHRES_ID,
        batch_key=BATCH_KEY,
    )
    sc_model = CondSCVI(sc_adata, weight_obs=False, prior="mog", num_classes_mog=10)
    sc_model.train(batch_size=1024, max_epochs=1)
    sc_model.save(VISIUM_COLON_SC_MODEL_PATH, overwrite=True, save_anndata=True)

    # Deconvolution with stLVM
    st_adata = st_adata[st_adata.layers["counts"].sum(1) > 10].copy()
    st_adata.obs[BATCH_KEY] = "spatial"

    def spatial_nn_gex_smth(stadata, n_neighs):
        sc.pp.neighbors(stadata, n_neighs, use_rep="spatial", key_added="Xspatial")
        stadata.obsp["Xspatial_connectivities"] = stadata.obsp["Xspatial_connectivities"].ceil()
        stadata.obsp["Xspatial_connectivities"].setdiag(1)
        return stadata.obsp["Xspatial_connectivities"].dot(stadata.layers["counts"])

    st_adata.layers["smoothed"] = spatial_nn_gex_smth(st_adata, n_neighs=5)
    st_model = DestVI.from_rna_model(
        st_adata,
        sc_model,
        add_celltypes=2,
        n_latent_amortization=None,
        anndata_setup_kwargs={"smoothed_layer": "smoothed"},
    )  # prior_mode = 'normal'
    st_model.train(
        max_epochs=1,
        n_epochs_kl_warmup=200,
        batch_size=1024,
        plan_kwargs={"weighting_kl_latent": 1e-2, "ct_sparsity_weight": 0},
    )
    st_model.save(VISIUM_COLON_ST_MODEL_PATH, overwrite=True)
    st_adata.obsm["proportions"] = st_model.get_proportions(keep_additional=True)
    ct_list = st_adata.obsm["proportions"].columns
    for ct in ct_list:
        data = st_adata.obsm["proportions"][ct].values
        st_adata.obs[ct] = data
    st_adata.obsm["fine_proportions"] = st_model.get_fine_celltypes(sc_model)
    ct_list = st_adata.obsm["fine_proportions"].columns
    for ct in ct_list:
        data = st_adata.obsm["fine_proportions"][ct].values
        st_adata.obs[ct] = data
    st_adata.var["betas"] = st_model.module.beta.detach().cpu().numpy()
    for ct, g in st_model.get_gamma().items():
        st_adata.obsm[f"{ct}_gamma"] = g
