import os

import numpy as np
import pandas as pd
import pytest

from scvi.data import synthetic_iid
from scvi.model import PEAKVI, SCANVI, SCVI, TOTALVI


def single_pass_for_online_update(model):
    dl = model._make_data_loader(model.adata, indices=range(0, 10))
    for tensors in dl:
        _, _, scvi_loss = model.module(tensors)
    scvi_loss.loss.backward()


def test_data_prep(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata2 has more genes and a perfect subset of adata1
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    SCVI.prepare_query_anndata(adata2, dir_path)
    SCVI.load_query_data(adata2, dir_path)

    adata3 = SCVI.prepare_query_anndata(adata2, dir_path, inplace=False)
    SCVI.load_query_data(adata3, dir_path)

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVI.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].X) == 0
    np.testing.assert_equal(
        adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy()
    )
    SCVI.load_query_data(adata4, dir_path)

    adata5 = SCVI.prepare_query_anndata(adata4, dir_path, inplace=False)
    SCVI.load_query_data(adata5, dir_path)


def test_data_prep_layer(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    adata1.layers["counts"] = adata1.X.copy()
    SCVI.setup_anndata(adata1, layer="counts", batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    adata4.layers["counts"] = adata4.X.copy()
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVI.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].layers["counts"]) == 0
    np.testing.assert_equal(
        adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy()
    )
    SCVI.load_query_data(adata4, dir_path)


def test_scvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, inplace_subset_query_vars=True)
    model2.train(max_epochs=1, plan_kwargs=dict(weight_decay=0.0))
    model2.get_latent_representation()

    # encoder linear layer equal
    one = (
        model.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    two = (
        model2.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    np.testing.assert_equal(one, two)
    single_pass_for_online_update(model2)
    assert (
        np.sum(
            model2.module.z_encoder.encoder.fc_layers[0][0]
            .weight.grad.cpu()
            .numpy()[:, : adata1.shape[1]]
        )
        == 0
    )
    # dispersion
    assert model2.module.px_r.requires_grad is False
    # library encoder linear layer
    assert model2.module.l_encoder.encoder.fc_layers[0][0].weight.requires_grad is True
    # 5 for n_latent, 4 for batches
    assert model2.module.decoder.px_decoder.fc_layers[0][0].weight.shape[1] == 9

    # test options
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(
        adata1,
        n_latent=n_latent,
        n_layers=2,
        encode_covariates=True,
        use_batch_norm="encoder",
        use_layer_norm="none",
    )
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(max_epochs=1, plan_kwargs=dict(weight_decay=0.0))
    # deactivate no grad decorator
    model2.get_latent_representation()
    # pytorch lightning zeros the grad, so this will get a grad to inspect
    single_pass_for_online_update(model2)
    grad = model2.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.sum(grad[:, -4:]) != 0
    grad = model2.module.decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) == 0

    # do not freeze expression
    model3 = SCVI.load_query_data(
        adata2,
        dir_path,
        freeze_expression=False,
        freeze_batchnorm_encoder=True,
        freeze_decoder_first_layer=False,
    )
    model3.train(max_epochs=1)
    model3.get_latent_representation()
    assert model3.module.z_encoder.encoder.fc_layers[0][1].momentum == 0
    # batch norm weight in encoder layer
    assert model3.module.z_encoder.encoder.fc_layers[0][1].weight.requires_grad is False
    single_pass_for_online_update(model3)
    grad = model3.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0
    grad = model3.module.decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0

    # do not freeze batchnorm
    model3 = SCVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=False)
    model3.train(max_epochs=1)
    model3.get_latent_representation()


def test_scvi_library_size_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent, use_observed_lib_size=False)

    assert (
        getattr(model.module, "library_log_means", None) is not None
        and model.module.library_log_means.shape == (1, 2)
        and model.module.library_log_means.count_nonzero().item() == 2
    )
    assert getattr(
        model.module, "library_log_vars", None
    ) is not None and model.module.library_log_vars.shape == (1, 2)

    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, inplace_subset_query_vars=True)
    assert (
        getattr(model2.module, "library_log_means", None) is not None
        and model2.module.library_log_means.shape == (1, 4)
        and model2.module.library_log_means[:, :2].equal(model.module.library_log_means)
        and model2.module.library_log_means.count_nonzero().item() == 4
    )
    assert (
        getattr(model2.module, "library_log_vars", None) is not None
        and model2.module.library_log_vars.shape == (1, 4)
        and model2.module.library_log_vars[:, :2].equal(model.module.library_log_vars)
    )


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

    # Test error on extra categoricals
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
    with pytest.raises(NotImplementedError):
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
    model2.train(max_epochs=1, plan_kwargs=dict(weight_decay=0.0))
    model2.get_latent_representation()
    model2.predict()

    # test classifier frozen
    class_query_weight = (
        model2.module.classifier.classifier[0]
        .fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()
    )
    class_ref_weight = (
        model.module.classifier.classifier[0]
        .fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()
    )
    # weight decay makes difference
    np.testing.assert_allclose(class_query_weight, class_ref_weight, atol=1e-07)

    # test classifier unfrozen
    model2 = SCANVI.load_query_data(adata2, dir_path, freeze_classifier=False)
    model2._unlabeled_indices = np.arange(adata2.n_obs)
    model2._labeled_indices = []
    model2.train(max_epochs=1)
    class_query_weight = (
        model2.module.classifier.classifier[0]
        .fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()
    )
    class_ref_weight = (
        model.module.classifier.classifier[0]
        .fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()
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


def test_totalvi_online_update(save_path):
    # basic case
    n_latent = 5
    adata1 = synthetic_iid()
    TOTALVI.setup_anndata(
        adata1,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    model = TOTALVI(adata1, n_latent=n_latent, use_batch_norm="decoder")
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = TOTALVI.load_query_data(adata2, dir_path)
    assert model2.module.background_pro_alpha.requires_grad is True
    model2.train(max_epochs=1)
    model2.get_latent_representation()

    # batch 3 has no proteins
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obsm["protein_expression"][adata2.obs.batch == "batch_3"] = 0

    # load from model in memory
    model3 = TOTALVI.load_query_data(adata2, model)
    model3.train(max_epochs=1)
    model3.get_latent_representation()


def test_peakvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    PEAKVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = PEAKVI(adata1, n_latent=n_latent)
    model.train(1, save_best=False)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = PEAKVI.load_query_data(adata2, dir_path)
    model2.train(max_epochs=1, weight_decay=0.0, save_best=False)
    model2.get_latent_representation()

    # encoder linear layer equal for peak features
    one = (
        model.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    two = (
        model2.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    np.testing.assert_equal(one, two)
    assert (
        np.sum(
            model2.module.z_encoder.encoder.fc_layers[0][0]
            .weight.grad.cpu()
            .numpy()[:, : adata1.shape[1]]
        )
        == 0
    )

    # test options
    adata1 = synthetic_iid()
    PEAKVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = PEAKVI(
        adata1,
        n_latent=n_latent,
        encode_covariates=True,
    )
    model.train(1, check_val_every_n_epoch=1, save_best=False)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = PEAKVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(max_epochs=1, weight_decay=0.0, save_best=False)
    # deactivate no grad decorator
    model2.get_latent_representation()
    # pytorch lightning zeros the grad, so this will get a grad to inspect
    single_pass_for_online_update(model2)
    grad = model2.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.count_nonzero(grad[:, -4:]) > 0

    # do not freeze expression
    model3 = PEAKVI.load_query_data(
        adata2,
        dir_path,
        freeze_expression=False,
        freeze_decoder_first_layer=False,
    )
    model3.train(max_epochs=1, save_best=False, weight_decay=0.0)
    model3.get_latent_representation()
    single_pass_for_online_update(model3)
    grad = model3.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.count_nonzero(grad[:, :-4]) != 0
    grad = model3.module.z_decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.count_nonzero(grad[:, :-4]) != 0
