import os

import numpy as np
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
    model.differential_expression(groupby="labels", group1="label_1", pseudocounts=3e-5)

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
    assert post_pred["rna"].shape == (n_obs, n_vars, 2)
    assert post_pred["protein"].shape == (n_obs, n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred["rna"].shape == (n_obs, n_vars)
    assert post_pred["protein"].shape == (n_obs, n_proteins)
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

    predictions, attributions = model.predict(ig_interpretability=True)
    # let's see an avg of score of top 5 markers for all samples put together
    ig_top_features = attributions.head(5)
    print(ig_top_features)

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
    assert totalanvi_pxr is not None
    assert totalvi_pxr is not None
    assert totalanvi_pxr is not totalvi_pxr
    totalanvi_model.train(1)

    totalanvi = TOTALANVI.from_totalvi_model(
        model, unlabeled_category="unlabeled", labels_key="labels", linear_classifier=True
    )
    totalanvi.train(
        max_epochs=1,
        adversarial_classifier=False,
        plan_kwargs={
            "pro_recons_weight": 0.3,
            "n_epochs_kl_warmup": 10.0,
            "lr": 3e-3,
            "classification_ratio": 1000.0,
            "max_kl_weight": 1.0,
        },
    )


def test_totalanvi_scarches(save_path):
    # test transfer_anndata_setup + view
    adata1 = synthetic_iid()
    TOTALANVI.setup_anndata(
        adata1,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        labels_key="labels",
        unlabeled_category="label_0",
    )
    model = TOTALANVI(adata1)
    model.train(1, train_size=0.5)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata2 has more genes and a perfect subset of adata1
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    TOTALANVI.prepare_query_anndata(adata2, dir_path)
    TOTALANVI_query = TOTALANVI.load_query_data(adata2, dir_path)
    TOTALANVI_query.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata3 = TOTALANVI.prepare_query_anndata(adata2, dir_path, inplace=False)
    TOTALANVI_query2 = TOTALANVI.load_query_data(adata3, dir_path)
    TOTALANVI_query2.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    TOTALANVI.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].X) == 0
    np.testing.assert_equal(adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy())
    TOTALANVI_query3 = TOTALANVI.load_query_data(adata4, dir_path)
    TOTALANVI_query3.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata5 = TOTALANVI.prepare_query_anndata(adata4, dir_path, inplace=False)
    TOTALANVI_query4 = TOTALANVI.load_query_data(adata5, dir_path)
    TOTALANVI_query4.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})
