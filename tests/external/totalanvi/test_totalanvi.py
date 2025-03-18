import pandas as pd
import pytest

from scvi.data import synthetic_iid
from scvi.external import TOTALANVI
from scvi.model import TOTALVI


def test_totalanvi():
    # test transfer_anndata_setup + view
    adata = synthetic_iid()
    TOTALANVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        labels_key="labels",
        unlabeled_category="label_0",
    )
    model = TOTALANVI(adata)
    model.train(1, train_size=0.5)
    model.differential_expression(groupby="labels", group1="label_1")

    adata = synthetic_iid()
    adata.obs["label_names"] = adata.obs["labels"]
    TOTALANVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        labels_key="labels",
        unlabeled_category="label_0",
    )

    n_obs = adata.n_obs
    n_vars = adata.n_vars
    n_proteins = adata.obsm["protein_expression"].shape[1]
    n_latent = 10

    model = TOTALANVI(adata, n_latent=n_latent)
    assert len(model._labeled_indices) == sum(adata.obs["labels"] != "label_0")
    assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == "label_0")
    model.train(1, train_size=0.5, check_val_every_n_epoch=1, batch_size=300)
    logged_keys = model.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "validation_classification_loss" in logged_keys
    assert "validation_accuracy" in logged_keys
    assert "validation_f1_score" in logged_keys
    assert "validation_calibration_error" in logged_keys

    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch=["batch_0", "batch_1"])
    model.get_latent_library_size()
    model.get_protein_foreground_probability()
    model.get_protein_foreground_probability(transform_batch=["batch_0", "batch_1"])
    post_pred = model.posterior_predictive_sample(n_samples=2)
    assert post_pred.shape == (n_obs, n_vars + n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred.shape == (n_obs, n_vars + n_proteins)
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(correlation_type="spearman")
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman", transform_batch=["batch_0", "batch_1"]
    )
    feature_correlation_matrix2 = model.get_feature_correlation_matrix(correlation_type="pearson")
    assert feature_correlation_matrix1.shape == (
        n_vars + n_proteins,
        n_vars + n_proteins,
    )
    assert feature_correlation_matrix2.shape == (
        n_vars + n_proteins,
        n_vars + n_proteins,
    )

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
    model.get_reconstruction_error(indices=model.validation_indices)

    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    df = model.predict(adata2, soft=True)
    assert isinstance(df, pd.DataFrame)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)

    # test from_totalvi_model
    adata = synthetic_iid()
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    model = TOTALVI(adata, n_latent=n_latent)
    a2 = synthetic_iid()
    totalanvi_model = TOTALANVI.from_totalvi_model(model, "label_0", labels_key="labels", adata=a2)
    with pytest.raises(ValueError):
        totalanvi_model = TOTALANVI.from_totalvi_model(model, "label_0", labels_key=None, adata=a2)

    # make sure the state_dicts are different objects for the two models
    assert totalanvi_model.module.state_dict() is not model.module.state_dict()
    totalanvi_pxr = totalanvi_model.module.state_dict().get("px_r", None)
    totalvi_pxr = model.module.state_dict().get("px_r", None)
    assert totalanvi_pxr is not None and totalvi_pxr is not None
    assert totalanvi_pxr is not totalvi_pxr
    totalanvi_model.train(1)
