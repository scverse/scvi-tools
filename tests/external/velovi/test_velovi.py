import scvelo as scv

from scvi.data import synthetic_iid
from scvi.external.velovi import VELOVI


def test_preprocess_data():
    adata = synthetic_iid()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    scv.pp.normalize_per_cell(adata)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    # TODO(adamgayoso): use real data for this test
    # preprocess_data(adata)


def test_velovi():
    n_latent = 5
    adata = synthetic_iid()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    VELOVI.setup_anndata(adata, unspliced_layer="unspliced", spliced_layer="spliced")
    model = VELOVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_velocity()
    model.get_latent_time()
    model.get_state_assignment()
    model.get_expression_fit()
    model.get_directional_uncertainty()
    model.get_permutation_scores(labels_key="labels")

    _ = model.history

    # tests __repr__
    print(model)
