import os

import numpy as np

from scvi.data import synthetic_iid
from scvi.model import AmortizedLDA


def test_lda_model_single_step(n_topics: int = 5):
    adata = synthetic_iid()
    AmortizedLDA.setup_anndata(adata)
    mod1 = AmortizedLDA(adata, n_topics=n_topics, cell_topic_prior=1.5, topic_feature_prior=1.5)
    mod1.train(max_steps=1, max_epochs=10)
    assert len(mod1.history["elbo_train"]) == 1


def test_lda_model_single_step_with_external_indices(n_topics: int = 5):
    adata = synthetic_iid()
    AmortizedLDA.setup_anndata(adata)
    mod1 = AmortizedLDA(adata, n_topics=n_topics, cell_topic_prior=1.5, topic_feature_prior=1.5)
    # in this case we will make a stratified version of indexing
    from sklearn.model_selection import train_test_split

    train_ind, valid_ind = train_test_split(
        adata.obs.batch.index.astype(int), test_size=0.6, stratify=adata.obs.batch
    )
    test_ind, valid_ind = train_test_split(
        valid_ind, test_size=0.5, stratify=adata.obs.batch[valid_ind]
    )
    mod1.train(
        max_steps=1,
        max_epochs=10,
        datasplitter_kwargs={
            "external_indexing": [np.array(train_ind), np.array(valid_ind), np.array(test_ind)]
        },
    )
    assert len(mod1.history["elbo_train"]) == 1


def test_lda_model(n_topics: int = 5):
    adata = synthetic_iid()

    # Test with float and Sequence priors.
    AmortizedLDA.setup_anndata(adata)
    mod1 = AmortizedLDA(adata, n_topics=n_topics, cell_topic_prior=1.5, topic_feature_prior=1.5)
    mod1.train(
        max_epochs=1,
        batch_size=256,
        lr=0.01,
    )
    mod2 = AmortizedLDA(
        adata,
        n_topics=n_topics,
        cell_topic_prior=[1.5 for _ in range(n_topics)],
        topic_feature_prior=[1.5 for _ in range(adata.n_vars)],
    )
    mod2.train(
        max_epochs=1,
        batch_size=256,
        lr=0.01,
    )

    mod = AmortizedLDA(adata, n_topics=n_topics)
    mod.train(
        max_epochs=5,
        batch_size=256,
        lr=0.01,
    )
    adata_gbt = mod.get_feature_by_topic().to_numpy()
    assert np.allclose(adata_gbt.sum(axis=0), 1)
    adata_lda = mod.get_latent_representation(adata).to_numpy()
    assert adata_lda.shape == (adata.n_obs, n_topics)
    assert np.all((adata_lda <= 1) & (adata_lda >= 0))
    assert np.allclose(adata_lda.sum(axis=1), 1)
    mod.get_elbo()
    mod.get_perplexity()

    adata2 = synthetic_iid()
    AmortizedLDA.setup_anndata(adata2)
    adata2_lda = mod.get_latent_representation(adata2).to_numpy()
    assert adata2_lda.shape == (adata2.n_obs, n_topics)
    assert np.all((adata2_lda <= 1) & (adata2_lda >= 0))
    assert np.allclose(adata2_lda.sum(axis=1), 1)
    mod.get_elbo(adata2)
    mod.get_perplexity(adata2)


def test_lda_model_save_load(save_path: str, n_topics: int = 5):
    adata = synthetic_iid()
    AmortizedLDA.setup_anndata(adata)
    mod = AmortizedLDA(adata, n_topics=n_topics)
    mod.train(
        max_epochs=5,
        batch_size=256,
        lr=0.01,
    )
    hist_elbo = mod.history_["elbo_train"]

    feature_by_topic_1 = mod.get_feature_by_topic(n_samples=5000)
    latent_1 = mod.get_latent_representation(n_samples=6000)

    save_path = os.path.join(save_path, "tmp")
    mod.save(save_path, overwrite=True, save_anndata=True)
    mod = AmortizedLDA.load(save_path)

    np.testing.assert_array_equal(mod.history_["elbo_train"], hist_elbo)

    feature_by_topic_2 = mod.get_feature_by_topic(n_samples=5000)
    latent_2 = mod.get_latent_representation(n_samples=6000)
    np.testing.assert_almost_equal(
        feature_by_topic_1.to_numpy(), feature_by_topic_2.to_numpy(), decimal=2
    )
    np.testing.assert_almost_equal(latent_1.to_numpy(), latent_2.to_numpy(), decimal=2)
