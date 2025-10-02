import numpy as np
import pandas as pd
import scipy
from scvi.data import synthetic_iid
from scvi.external import SCWAE


def test_scwae():
    adata = synthetic_iid()
    adata.X = scipy.sparse.csr_matrix(adata.X)
    SCWAE.setup_anndata(adata)

    model = SCWAE(adata)
    model.train(
        1, check_val_every_n_epoch=1, train_size=0.5,
        plan_kwargs={'wae_lambda': 1., 'weighting_random_encoder': 1.})
    model.get_latent_representation()

    # tests __repr__
    print(model)


def test_scwae_scarches():
    # reference adata
    adata = synthetic_iid()
    SCWAE.setup_anndata(adata)
    model = SCWAE(adata)
    model.train(max_epochs=1)

    # Query adata
    adata2 = synthetic_iid()

    # Make query adata and model
    SCWAE.prepare_query_anndata(adata2, model)
    model2 = SCWAE.load_query_data(adata2, model)
    model2.train(max_epochs=1)
    assert (
        model2.get_latent_representation(
            adata=adata,
        ).shape[0]
        == adata.shape[0]
    )
    assert (
        model2.get_latent_representation(
            adata=adata2,
        ).shape[0]
        == adata2.shape[0]
    )
