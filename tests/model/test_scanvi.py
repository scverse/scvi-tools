import os
import pickle

import numpy as np
import pandas as pd
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.model import SCANVI, SCVI
from scvi.utils import attrdict


def test_saving_and_loading(save_path):
    def legacy_save(
        model,
        dir_path,
        prefix=None,
        overwrite=False,
        save_anndata=False,
        **anndata_write_kwargs,
    ):
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide an unexisting directory for saving."
            )

        file_name_prefix = prefix or ""

        if save_anndata:
            model.adata.write(
                os.path.join(dir_path, f"{file_name_prefix}adata.h5ad"),
                **anndata_write_kwargs,
            )

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
        attr_save_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
        varnames_save_path = os.path.join(dir_path, f"{file_name_prefix}var_names.csv")

        torch.save(model.module.state_dict(), model_save_path)

        var_names = model.adata.var_names.astype(str)
        var_names = var_names.to_numpy()
        np.savetxt(varnames_save_path, var_names, fmt="%s")

        # get all the user attributes
        user_attributes = model._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        with open(attr_save_path, "wb") as f:
            pickle.dump(user_attributes, f)

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

    prefix = "SCANVI_"
    SCANVI.setup_anndata(adata, "labels", "label_0", batch_key="batch")
    model = SCANVI(adata)
    model.train(max_epochs=1, train_size=0.5)
    p1 = model.predict()
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model.view_setup_args(save_path, prefix=prefix)
    model = SCANVI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        SCANVI.load(save_path, adata=tmp_adata, prefix=prefix)
    model = SCANVI.load(save_path, adata=adata, prefix=prefix)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # Test legacy loading
    legacy_save_path = os.path.join(save_path, "legacy/")
    legacy_save(model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix)
    with pytest.raises(ValueError):
        SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    SCANVI.convert_legacy_save(legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix)
    m = SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)


def test_scanvi():
    adata = synthetic_iid()
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    model = SCANVI(adata, n_latent=10)
    assert len(model._labeled_indices) == sum(adata.obs["labels"] != "label_0")
    assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == "label_0")
    model.train(1, train_size=0.5, check_val_every_n_epoch=1)
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
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    df = model.predict(adata2, soft=True)
    assert isinstance(df, pd.DataFrame)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")

    # test that all data labeled runs
    unknown_label = "asdf"
    a = synthetic_iid()
    SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = SCANVI(a)
    m.train(1)

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    a = synthetic_iid()
    SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = SCANVI(a)
    m.train(1, train_size=0.9)

    # test from_scvi_model
    a = synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
    )
    m = SCVI(a, use_observed_lib_size=False)
    a2 = synthetic_iid()
    scanvi_model = SCANVI.from_scvi_model(m, "label_0", labels_key="labels", adata=a2)
    with pytest.raises(ValueError):
        scanvi_model = SCANVI.from_scvi_model(m, "label_0", labels_key=None, adata=a2)

    # make sure the state_dicts are different objects for the two models
    assert scanvi_model.module.state_dict() is not m.module.state_dict()
    scanvi_pxr = scanvi_model.module.state_dict().get("px_r", None)
    scvi_pxr = m.module.state_dict().get("px_r", None)
    assert scanvi_pxr is not None
    assert scvi_pxr is not None
    assert scanvi_pxr is not scvi_pxr
    scanvi_model.train(1)

    # Test without label groups
    scanvi_model = SCANVI.from_scvi_model(
        m, "label_0", labels_key="labels", use_labels_groups=False
    )
    scanvi_model.train(1)

    # test from_scvi_model with size_factor
    a = synthetic_iid()
    a.obs["size_factor"] = np.random.randint(1, 5, size=(a.shape[0],))
    SCVI.setup_anndata(a, batch_key="batch", labels_key="labels", size_factor_key="size_factor")
    m = SCVI(a, use_observed_lib_size=False)
    a2 = synthetic_iid()
    a2.obs["size_factor"] = np.random.randint(1, 5, size=(a2.shape[0],))
    scanvi_model = SCANVI.from_scvi_model(m, "label_0", adata=a2)
    scanvi_model.train(1)


def test_scanvi_with_external_indices():
    adata = synthetic_iid()
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    model = SCANVI(adata, n_latent=10)
    assert len(model._labeled_indices) == sum(adata.obs["labels"] != "label_0")
    assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == "label_0")
    # in this case we will make a stratified version of indexing
    from sklearn.model_selection import train_test_split

    train_ind, valid_ind = train_test_split(
        adata.obs.batch.index.astype(int), test_size=0.6, stratify=adata.obs.batch
    )
    test_ind, valid_ind = train_test_split(
        valid_ind, test_size=0.5, stratify=adata.obs.batch[valid_ind]
    )
    model.train(
        1,
        check_val_every_n_epoch=1,
        datasplitter_kwargs={
            "external_indexing": [np.array(train_ind), np.array(valid_ind), np.array(test_ind)]
        },
    )
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
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    df = model.predict(adata2, soft=True)
    assert isinstance(df, pd.DataFrame)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")

    # test that all data labeled runs
    unknown_label = "asdf"
    a = synthetic_iid()
    SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = SCANVI(a)
    m.train(1)

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    a = synthetic_iid()
    SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = SCANVI(a)
    m.train(1, train_size=0.9)

    # test from_scvi_model
    a = synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
    )
    m = SCVI(a, use_observed_lib_size=False)
    a2 = synthetic_iid()
    scanvi_model = SCANVI.from_scvi_model(m, "label_0", labels_key="labels", adata=a2)
    with pytest.raises(ValueError):
        scanvi_model = SCANVI.from_scvi_model(m, "label_0", labels_key=None, adata=a2)

    # make sure the state_dicts are different objects for the two models
    assert scanvi_model.module.state_dict() is not m.module.state_dict()
    scanvi_pxr = scanvi_model.module.state_dict().get("px_r", None)
    scvi_pxr = m.module.state_dict().get("px_r", None)
    assert scanvi_pxr is not None
    assert scvi_pxr is not None
    assert scanvi_pxr is not scvi_pxr
    scanvi_model.train(1)

    # Test without label groups
    scanvi_model = SCANVI.from_scvi_model(
        m, "label_0", labels_key="labels", use_labels_groups=False
    )
    scanvi_model.train(1)

    # test from_scvi_model with size_factor
    a = synthetic_iid()
    a.obs["size_factor"] = np.random.randint(1, 5, size=(a.shape[0],))
    SCVI.setup_anndata(a, batch_key="batch", labels_key="labels", size_factor_key="size_factor")
    m = SCVI(a, use_observed_lib_size=False)
    a2 = synthetic_iid()
    a2.obs["size_factor"] = np.random.randint(1, 5, size=(a2.shape[0],))
    scanvi_model = SCANVI.from_scvi_model(m, "label_0", adata=a2)
    scanvi_model.train(1)


def test_scanvi_predict_use_posterior_mean():
    adata = synthetic_iid()
    SCANVI.setup_anndata(adata, labels_key="labels", unlabeled_category="label_0")

    model = SCANVI(adata)
    model.train(max_epochs=1)

    _ = model.predict(use_posterior_mean=True)
    _ = model.predict(use_posterior_mean=False)


def test_linear_classifier_scanvi(n_latent: int = 10, n_labels: int = 5):
    adata = synthetic_iid(n_labels=n_labels)
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    model = SCANVI(adata, linear_classifier=True, n_latent=n_latent)

    assert isinstance(model.module.classifier.classifier[0], torch.nn.Linear)
    assert model.module.classifier.classifier[0].in_features == n_latent
    assert model.module.classifier.classifier[0].out_features == n_labels - 1

    model.train(max_epochs=1)


def test_multiple_covariates_scanvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

    SCANVI.setup_anndata(
        adata,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCANVI(adata)
    m.train(1)
    m.get_latent_representation()
    m.get_elbo()
    m.get_marginal_ll(n_mc_samples=3)
    m.get_marginal_ll(adata, return_mean=True, n_mc_samples=6, n_mc_samples_per_pass=1)
    m.get_marginal_ll(adata, return_mean=True, n_mc_samples=6, n_mc_samples_per_pass=6)
    m.differential_expression(
        idx1=np.arange(50), idx2=51 + np.arange(50), mode="vanilla", weights="uniform"
    )
    m.differential_expression(
        idx1=np.arange(50),
        idx2=51 + np.arange(50),
        mode="vanilla",
        weights="importance",
        importance_weighting_kwargs={"n_mc_samples": 10, "n_mc_samples_per_pass": 1},
    )
    m.differential_expression(
        idx1=np.arange(50),
        idx2=51 + np.arange(50),
        mode="vanilla",
        weights="importance",
        importance_weighting_kwargs={"n_mc_samples": 10, "n_mc_samples_per_pass": 10},
    )
    m.get_reconstruction_error()
    m.get_normalized_expression(n_samples=1)
    m.get_normalized_expression(n_samples=2)


def test_multiple_encoded_covariates_scanvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

    SCANVI.setup_anndata(
        adata,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCANVI(adata, encode_covariates=True)
    m.train(1)


def test_scanvi_online_update(save_path):
    # ref has semi-observed labels
    n_latent = 5
    adata1 = synthetic_iid()
    new_labels = adata1.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata1.obs["labels"] = pd.Categorical(new_labels)
    adata1.obs["cont1"] = np.random.normal(size=(adata1.shape[0],))
    adata1.obs["cont2"] = np.random.normal(size=(adata1.shape[0],))

    SCANVI.setup_anndata(
        adata1,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
    )
    model = SCANVI(
        adata1,
        n_latent=n_latent,
        encode_covariates=True,
    )
    model.train(max_epochs=1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # query has all missing labels
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obs["labels"] = "Unknown"
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))

    model = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.predict()

    # query has all missing labels and no labels key
    del adata2.obs["labels"]

    model = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.predict()

    # query has no missing labels
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))

    model = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.predict()

    # Test on extra categoricals as well
    adata1 = synthetic_iid()
    new_labels = adata1.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata1.obs["labels"] = pd.Categorical(new_labels)
    adata1.obs["cont1"] = np.random.normal(size=(adata1.shape[0],))
    adata1.obs["cont2"] = np.random.normal(size=(adata1.shape[0],))
    adata1.obs["cat1"] = np.random.randint(0, 5, size=(adata1.shape[0],))
    adata1.obs["cat2"] = np.random.randint(0, 5, size=(adata1.shape[0],))
    SCANVI.setup_anndata(
        adata1,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    model = SCANVI(
        adata1,
        n_latent=n_latent,
        encode_covariates=True,
    )
    model.train(max_epochs=1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obs["labels"] = "Unknown"
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cat1"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    adata2.obs["cat2"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)

    # ref has fully-observed labels
    n_latent = 5
    adata1 = synthetic_iid()
    new_labels = adata1.obs.labels.to_numpy()
    adata1.obs["labels"] = pd.Categorical(new_labels)
    SCANVI.setup_anndata(adata1, "labels", "Unknown", batch_key="batch")
    model = SCANVI(adata1, n_latent=n_latent, encode_covariates=True)
    model.train(max_epochs=1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # query has one new label
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    new_labels = adata2.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata2.obs["labels"] = pd.Categorical(new_labels)

    model2 = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model2._unlabeled_indices = np.arange(adata2.n_obs)
    model2._labeled_indices = []
    model2.train(max_epochs=1, plan_kwargs={"weight_decay": 0.0})
    model2.get_latent_representation()
    model2.predict()

    # test classifier frozen
    class_query_weight = (
        model2.module.classifier.classifier[0].fc_layers[0][0].weight.detach().cpu().numpy()
    )
    class_ref_weight = (
        model.module.classifier.classifier[0].fc_layers[0][0].weight.detach().cpu().numpy()
    )
    # weight decay makes difference
    np.testing.assert_allclose(class_query_weight, class_ref_weight, atol=1e-07)

    # test classifier unfrozen
    model2 = SCANVI.load_query_data(adata2, dir_path, freeze_classifier=False)
    model2._unlabeled_indices = np.arange(adata2.n_obs)
    model2._labeled_indices = []
    model2.train(max_epochs=1)
    class_query_weight = (
        model2.module.classifier.classifier[0].fc_layers[0][0].weight.detach().cpu().numpy()
    )
    class_ref_weight = (
        model.module.classifier.classifier[0].fc_layers[0][0].weight.detach().cpu().numpy()
    )
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(class_query_weight, class_ref_weight, atol=1e-07)

    # test saving and loading of online scanvi
    a = synthetic_iid()
    ref = a[a.obs["labels"] != "label_2"].copy()  # only has labels 0 and 1
    SCANVI.setup_anndata(ref, "labels", "label_2", batch_key="batch")
    m = SCANVI(ref)
    m.train(max_epochs=1)
    m.save(save_path, overwrite=True)
    query = a[a.obs["labels"] != "label_0"].copy()
    query = synthetic_iid()  # has labels 0 and 2. 2 is unknown
    m_q = SCANVI.load_query_data(query, save_path)
    m_q.save(save_path, overwrite=True)
    m_q = SCANVI.load(save_path, adata=query)
    m_q.predict()
    m_q.get_elbo()


def test_scanvi_logits_backwards_compat(save_path: str):
    """Check that we can replicate pre-fix behavior with `logits=False`."""
    adata = synthetic_iid()
    SCANVI.setup_anndata(adata, labels_key="labels", unlabeled_category="label_0")

    model = SCANVI(adata, classifier_parameters={"logits": False})
    model.train(max_epochs=1)

    model_path = os.path.join(save_path, "scanvi_logits_model")
    model.save(model_path, overwrite=True)
    del model

    model = SCANVI.load(model_path, adata)
    assert not model.module.classifier.logits
    assert isinstance(model.module.classifier.classifier[-1], torch.nn.Softmax)


def test_scanvi_pre_logits_fix_load(save_path: str):
    """See #2310. Check old model saves use the old behavior."""
    resave_model_path = os.path.join(save_path, "pre_logits_fix_scanvi")
    model = SCANVI.load(resave_model_path)

    def check_no_logits_and_softmax(model: SCANVI):
        assert not model.module.classifier.logits
        assert isinstance(model.module.classifier.classifier[-1], torch.nn.Softmax)

    check_no_logits_and_softmax(model)

    model.save(resave_model_path, overwrite=True)
    adata = model.adata
    del model

    model = SCANVI.load(resave_model_path, adata)
    check_no_logits_and_softmax(model)


@pytest.mark.parametrize("unlabeled_cat", ["label_0"])
def test_scanvi_interpretability_ig(unlabeled_cat: str):
    adata = synthetic_iid(batch_size=50)
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    SCANVI.setup_anndata(
        adata,
        labels_key="labels",
        unlabeled_category=unlabeled_cat,
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    model = SCANVI(adata, n_latent=10)
    model.train(1, train_size=0.5, check_val_every_n_epoch=1)

    # get the IG for all data
    predictions, attributions = model.predict(ig_interpretability=True)  # orignal predictions
    # let's see an avg of score of top 5 genes for all samples put together
    ig_top_features = attributions.head(5)
    print(ig_top_features)

    # new data ig prediction specific for samples, top 5 genes
    adata2 = synthetic_iid(batch_size=10)
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cat1"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    adata2.obs["cat2"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    predictions, attributions = model.predict(adata2, indices=[1, 2, 3], ig_interpretability=True)
    ig_top_features_3_samples = attributions.head(5)
    print(ig_top_features_3_samples)


@pytest.mark.parametrize("unlabeled_cat", ["label_0"])
def test_scanvi_interpretability_shap(unlabeled_cat: str):
    adata = synthetic_iid(batch_size=50)
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    SCANVI.setup_anndata(
        adata,
        labels_key="labels",
        unlabeled_category=unlabeled_cat,
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    model = SCANVI(adata, n_latent=10)
    model.train(1, train_size=0.5, check_val_every_n_epoch=1)

    # new data for shap prediction specific for samples, top 5 genes
    adata2 = synthetic_iid(batch_size=10)
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cat1"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    adata2.obs["cat2"] = np.random.randint(0, 5, size=(adata2.shape[0],))

    # now run shap values and compare to previous results
    # (here, the more labels the more time it will take to run)
    shap_values = model.shap_predict(shap_args={"nsamples": 100})
    # select the label we want to understand (usually the '1' class)
    shap_top_features = model.get_ranked_markers(attrs=shap_values[:, :, 1]).head(5)
    print(shap_top_features)

    # now run shap values for the test set (can be specific class or indices and with params)
    # (here, the more labels the more time it will take to run)
    shap_values_test = model.shap_predict(
        adata2,
        indices=[1, 2, 3],
        shap_args={"link": "identity", "silent": True, "gc_collect": True, "nsamples": 300},
    )
    # # select the label we want to understand (usually the '1' class)
    shap_top_features_test = model.get_ranked_markers(attrs=shap_values_test[:, :, 1]).head(5)
    print(shap_top_features_test)
