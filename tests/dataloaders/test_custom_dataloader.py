from __future__ import annotations

import os
from pprint import pprint

import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.dataloaders import MappedCollectionDataModule, TileDBDataModule
from scvi.external import MRVI
from scvi.utils import dependencies


@pytest.fixture(scope="module")
def setup_lamindb_instance():
    os.system("lamin init --storage ./test_lamindb_loader")
    yield
    os.system("rm -r ./test_lamindb_loader")
    os.system("lamin delete --force test_lamindb_loader")


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_small(save_path: str, setup_lamindb_instance):
    import lamindb as ln

    # prepare test data
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()

    artifact1 = ln.Artifact.from_anndata(adata1, key="part_one.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="part_two.h5ad").save()

    collection = ln.Collection([artifact1, artifact2], key="gather")
    collection.save()

    artifacts = collection.artifacts.all()
    artifacts.df()

    datamodule = MappedCollectionDataModule(
        collection,
        batch_key="batch",
        batch_size=1024,
        join="inner",
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    # pprint(datamodule.registry)

    model = scvi.model.SCVI(registry=datamodule.registry)
    # pprint(model.summary_stats)
    # pprint(model.module)

    model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    model.history.keys()

    # The way to extract the internal model analysis is by the inference_dataloader
    # Datamodule will always require passing it into all downstream functions.
    inference_dataloader = datamodule.inference_dataloader()
    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)
    _ = model.posterior_predictive_sample(dataloader=inference_dataloader)
    _ = model.get_normalized_expression(dataloader=inference_dataloader)
    _ = model.get_likelihood_parameters(dataloader=inference_dataloader)
    _ = model._get_denoised_samples(dataloader=inference_dataloader)
    _ = model.get_latent_library_size(dataloader=inference_dataloader, give_mean=False)

    # repeat but with other data with fewer indices and smaller batch size
    adata1_small = synthetic_iid(batch_size=10)
    adata2_small = synthetic_iid(batch_size=10)
    artifact1_small = ln.Artifact.from_anndata(adata1_small, key="part_one_small.h5ad").save()
    artifact2_small = ln.Artifact.from_anndata(adata2_small, key="part_two_small.h5ad").save()
    collection_small = ln.Collection([artifact1_small, artifact2_small], key="gather")
    collection_small.save()
    datamodule_small = MappedCollectionDataModule(
        collection_small,
        batch_key="batch",
        batch_size=1024,
        join="inner",
        collection_val=collection,
    )
    inference_dataloader_small = datamodule_small.inference_dataloader(batch_size=128)
    _ = model.get_elbo(return_mean=False, dataloader=inference_dataloader_small)
    _ = model.get_marginal_ll(n_mc_samples=3, dataloader=inference_dataloader_small)
    _ = model.get_reconstruction_error(return_mean=False, dataloader=inference_dataloader_small)
    _ = model.get_latent_representation(dataloader=inference_dataloader_small)
    _ = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["gene_1", "gene_2"], dataloader=inference_dataloader_small
    )
    _ = model.get_normalized_expression(n_samples=2, dataloader=inference_dataloader_small)

    # load and save and make a query with the other data
    model.save("lamin_model", save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis
    loaded_model = scvi.model.SCVI.load("lamin_model", adata=False)
    loaded_model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    model_query = model.load_query_data(
        adata=False, reference_model="lamin_model", registry=datamodule.registry
    )
    model_query.train(
        max_epochs=1, datamodule=datamodule_small, check_val_every_n_epoch=1, train_size=0.9
    )
    model_query.history.keys()

    _ = model_query.get_elbo(dataloader=inference_dataloader_small)
    _ = model_query.get_marginal_ll(dataloader=inference_dataloader_small)
    _ = model_query.get_reconstruction_error(dataloader=inference_dataloader_small)
    _ = model_query.get_latent_representation(dataloader=inference_dataloader_small)
    _ = model_query.posterior_predictive_sample(dataloader=inference_dataloader_small)
    _ = model_query.get_normalized_expression(dataloader=inference_dataloader_small)
    _ = model_query.get_likelihood_parameters(dataloader=inference_dataloader_small)
    _ = model_query._get_denoised_samples(dataloader=inference_dataloader_small)
    _ = model_query.get_latent_library_size(dataloader=inference_dataloader_small, give_mean=False)

    # query again but with the adata of the model, which might bring more functionality
    adata = collection.load(join="inner")
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    with pytest.raises(ValueError):
        model.load_query_data(adata=adata)
    model_query_adata = model.load_query_data(adata=adata, reference_model="lamin_model")
    model_query_adata.train(max_epochs=1, check_val_every_n_epoch=1, train_size=0.9)
    model_query_adata.history.keys()
    _ = model_query_adata.get_elbo()
    _ = model_query_adata.get_marginal_ll()
    _ = model_query_adata.get_reconstruction_error()
    _ = model_query_adata.get_latent_representation()
    _ = model_query_adata.get_latent_representation(dataloader=inference_dataloader)
    _ = model_query_adata.posterior_predictive_sample(indices=[1, 2, 3])
    model.save("lamin_model", save_anndata=False, overwrite=True, datamodule=datamodule)
    model.load("lamin_model", adata=False)
    model.load_query_data(adata=False, reference_model="lamin_model", registry=datamodule.registry)

    # create a regular model
    model.load_query_data(adata=adata, reference_model="lamin_model")
    model_adata = model.load("lamin_model", adata=adata)
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model_adata.train(max_epochs=1, check_val_every_n_epoch=1, train_size=0.9)
    model_adata.save(
        "lamin_model_anndata", save_anndata=True, overwrite=True, datamodule=datamodule
    )
    model_adata.load("lamin_model_anndata")
    model_adata.load_query_data(
        adata=adata, reference_model="lamin_model_anndata", registry=datamodule.registry
    )
    model_adata.history.keys()
    # test different gene_likelihoods
    for gene_likelihood in ["zinb", "nb", "poisson"]:
        model_adata = scvi.model.SCVI(adata, gene_likelihood=gene_likelihood)
        model_adata.train(1, check_val_every_n_epoch=1, train_size=0.9)
        model_adata.posterior_predictive_sample()
        model_adata.get_latent_representation()
        model_adata.get_normalized_expression()

    scvi.model.SCVI.prepare_query_anndata(adata, "lamin_model", return_reference_var_names=True)
    scvi.model.SCVI.load_query_data(adata, "lamin_model")

    scvi.model.SCVI.prepare_query_anndata(adata, model_query)
    scvi.model.SCVI.load_query_data(adata, model_query_adata)


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scanvi_small(save_path: str, setup_lamindb_instance):
    import lamindb as ln

    # prepare test data
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()

    artifact1_scanvi = ln.Artifact.from_anndata(adata1, key="part_one_scanvi.h5ad").save()
    artifact2_scanvi = ln.Artifact.from_anndata(adata2, key="part_two_scanvi.h5ad").save()

    collection_scanvi = ln.Collection([artifact1_scanvi, artifact2_scanvi], key="gather")
    collection_scanvi.save()

    artifacts = collection_scanvi.artifacts.all()
    artifacts.df()

    datamodule = MappedCollectionDataModule(
        collection_scanvi,
        label_key="labels",
        batch_key="batch",
        batch_size=1024,
        join="inner",
        model_name="SCANVI",
    )

    # We can now create the scVI model object and train it:
    model = scvi.model.SCANVI(
        adata=None,
        registry=datamodule.registry,
        encode_covariates=False,
        datamodule=datamodule,
    )

    model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=1,
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)
    model.history.keys()

    # save the model
    model.save(
        "lamin_model_scanvi",
        save_anndata=False,
        overwrite=True,
        datamodule=datamodule,
    )
    # load it back and do downstream analysis
    loaded_model = scvi.model.SCANVI.load("lamin_model_scanvi", adata=False, datamodule=datamodule)
    loaded_model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=1,
        early_stopping=False,
    )

    inference_dataloader = datamodule.inference_dataloader()

    model.predict(dataloader=inference_dataloader, soft=False)

    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)
    _ = model.posterior_predictive_sample(dataloader=inference_dataloader)
    _ = model.get_normalized_expression(dataloader=inference_dataloader)
    _ = model.get_likelihood_parameters(dataloader=inference_dataloader)
    _ = model._get_denoised_samples(dataloader=inference_dataloader)
    _ = model.get_latent_library_size(dataloader=inference_dataloader, give_mean=False)
    _ = loaded_model.get_elbo(dataloader=inference_dataloader)

    logged_keys = model.history.keys()
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "train_classification_loss" in logged_keys
    assert "train_accuracy" in logged_keys
    assert "train_f1_score" in logged_keys
    assert "train_calibration_error" in logged_keys

    # repeat but with other data with fewer indices and smaller batch size
    adata1_scanvi_small = synthetic_iid(batch_size=10)
    adata2_scanvi_small = synthetic_iid(batch_size=10)
    artifact1_scanvi_small = ln.Artifact.from_anndata(
        adata1_scanvi_small, key="part_one_small_scanvi.h5ad"
    ).save()
    artifact2_scanvi_small = ln.Artifact.from_anndata(
        adata2_scanvi_small, key="part_two_small_scanvi.h5ad"
    ).save()
    collection_scanvi_small = ln.Collection(
        [artifact1_scanvi_small, artifact2_scanvi_small], key="gather"
    )
    collection_scanvi_small.save()
    datamodule_small = MappedCollectionDataModule(
        collection_scanvi_small,
        batch_key="batch",
        batch_size=1024,
        join="inner",
        model_name="SCANVI",
    )
    inference_dataloader_small = datamodule_small.inference_dataloader(batch_size=128)

    model.predict(dataloader=inference_dataloader_small, soft=False)

    # train from scvi model
    model_scvi = scvi.model.SCVI(registry=datamodule.registry)

    # with validation collection
    datamodule = MappedCollectionDataModule(
        collection_scanvi,
        label_key="labels",
        batch_key="batch",
        batch_size=1024,
        join="inner",
        collection_val=collection_scanvi_small,
    )

    model_scvi.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
        check_val_every_n_epoch=1,
        train_size=0.9,
    )
    model_scvi.save("lamin_model_scvi", save_anndata=False, overwrite=True, datamodule=datamodule)

    logged_keys = model_scvi.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "validation_loss" in logged_keys

    # We can now create the scVI model object and train it:
    model_scanvi_from_scvi = scvi.model.SCANVI.from_scvi_model(
        scvi_model=model_scvi,
        adata=None,
        registry=datamodule.registry,
        encode_covariates=False,
        unlabeled_category="label_0",
        labels_key="labels",
        datamodule=datamodule,
    )
    model_scanvi_from_scvi.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )
    # save the model
    model_scanvi_from_scvi.save(
        "lamin_model_scanvi_from_scvi", save_anndata=False, overwrite=True, datamodule=datamodule
    )
    # load it back and do downstream analysis
    scvi.model.SCANVI.load("lamin_model_scanvi_from_scvi", adata=False, datamodule=datamodule)

    logged_keys = model_scanvi_from_scvi.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "validation_loss" in logged_keys

    inference_dataloader = datamodule.inference_dataloader()

    _ = model_scanvi_from_scvi.get_elbo(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_marginal_ll(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_latent_representation(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.posterior_predictive_sample(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_normalized_expression(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_likelihood_parameters(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi._get_denoised_samples(dataloader=inference_dataloader)
    _ = model_scanvi_from_scvi.get_latent_library_size(
        dataloader=inference_dataloader, give_mean=False
    )

    # create scanvi from adata
    adata = collection_scanvi.load(join="inner")  # we can continue to

    scvi.model.SCANVI.setup_anndata(
        adata, batch_key="batch", labels_key="labels", unlabeled_category="label_0"
    )
    model_query_adata = scvi.model.SCANVI(adata, encode_covariates=True)
    model_query_adata.train(max_epochs=1, check_val_every_n_epoch=1, train_size=0.9)
    model_query_adata.predict(adata=adata, soft=True)
    model.history.keys()


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_mrvi_small(save_path: str, setup_lamindb_instance):
    import lamindb as ln

    # prepare test data
    adata1 = synthetic_iid()
    adata1.obs.index.name = "cell_id"
    adata1.obs["sample"] = np.random.choice(15, size=adata1.shape[0])
    adata1.obs["sample_str"] = [chr(i + ord("a")) for i in adata1.obs["sample"]]
    adata1.obs["sample_str"] = adata1.obs["sample_str"].astype(str)
    meta11 = np.random.randint(0, 2, size=15)
    adata1.obs["meta1"] = meta11[adata1.obs["sample"].values]
    meta21 = np.random.randn(15)
    adata1.obs["meta2"] = meta21[adata1.obs["sample"].values]
    adata1.obs["cont_cov"] = np.random.normal(0, 1, size=adata1.shape[0])
    adata1.obs["meta1_cat"] = "CAT_" + adata1.obs["meta1"].astype(str)
    adata1.obs.loc[:, "disjoint_batch"] = (adata1.obs.loc[:, "sample"] <= 6).replace(
        {True: "batch_0", False: "batch_1"}
    )
    adata1.obs["dummy_batch"] = 1

    adata2 = synthetic_iid()
    adata2.obs.index.name = "cell_id"
    adata2.obs["sample"] = np.random.choice(15, size=adata2.shape[0])
    adata2.obs["sample_str"] = [chr(i + ord("a")) for i in adata2.obs["sample"]]
    adata2.obs["sample_str"] = adata2.obs["sample_str"].astype(str)
    meta12 = np.random.randint(0, 2, size=15)
    adata2.obs["meta1"] = meta12[adata2.obs["sample"].values]
    meta22 = np.random.randn(15)
    adata2.obs["meta2"] = meta22[adata2.obs["sample"].values]
    adata2.obs["cont_cov"] = np.random.normal(0, 1, size=adata2.shape[0])
    adata2.obs["meta1_cat"] = "CAT_" + adata2.obs["meta1"].astype(str)
    adata2.obs.loc[:, "disjoint_batch"] = (adata2.obs.loc[:, "sample"] <= 6).replace(
        {True: "batch_0", False: "batch_1"}
    )
    adata2.obs["dummy_batch"] = 2

    artifact1 = ln.Artifact.from_anndata(adata1, key="part_one.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="part_two.h5ad").save()

    collection = ln.Collection([artifact1, artifact2], key="gather")
    collection.save()

    artifacts = collection.artifacts.all()
    artifacts.df()

    datamodule = MappedCollectionDataModule(
        collection,
        batch_key="batch",
        sample_key="sample_str",
        batch_size=1024,
        join="inner",
        model_name="TorchMRVI",
        collection_val=collection,
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    # pprint(datamodule.registry)

    model = MRVI(registry=datamodule.registry)
    # pprint(model.summary_stats)
    # pprint(model.module)

    model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    logged_keys = model.history.keys()
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "kl_global_train" in logged_keys
    assert "validation_loss" in logged_keys
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "kl_global_validation" in logged_keys

    # The way to extract the internal model analysis is by the inference_dataloader
    # Datamodule will always require passing it into all downstream functions.
    inference_dataloader = datamodule.inference_dataloader()
    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(give_z=False, dataloader=inference_dataloader)
    _ = model.get_latent_representation(give_z=True, dataloader=inference_dataloader)

    # save and load model
    model.save(
        "lamin_model_mrvi",
        save_anndata=False,
        overwrite=True,
        datamodule=datamodule,
    )
    # load it back and do downstream analysis
    loaded_model = MRVI.load("lamin_model_mrvi", adata=False)
    loaded_model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=1,
        early_stopping=False,
    )

    _ = loaded_model.get_elbo(dataloader=inference_dataloader)
    _ = loaded_model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = loaded_model.get_latent_representation(give_z=False, dataloader=inference_dataloader)
    _ = loaded_model.get_latent_representation(give_z=True, dataloader=inference_dataloader)

    loaded_logged_keys = loaded_model.history.keys()
    assert "elbo_train" in loaded_logged_keys
    assert "reconstruction_loss_train" in loaded_logged_keys
    assert "kl_local_train" in loaded_logged_keys
    assert "kl_global_train" in loaded_logged_keys
    assert "validation_loss" in loaded_logged_keys
    assert "elbo_validation" in loaded_logged_keys
    assert "reconstruction_loss_validation" in loaded_logged_keys
    assert "kl_local_validation" in loaded_logged_keys
    assert "kl_global_validation" in loaded_logged_keys


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_small_with_covariates(save_path: str, setup_lamindb_instance):
    import lamindb as ln

    # prepare test data
    adata1 = synthetic_iid()
    adata1.obs["cat1"] = np.random.randint(0, 5, size=(adata1.shape[0],))
    adata1.obs["cat2"] = np.random.randint(0, 5, size=(adata1.shape[0],))
    adata1.obs["cont1"] = np.random.normal(size=(adata1.shape[0],))
    adata1.obs["cont2"] = np.random.normal(size=(adata1.shape[0],))
    adata2 = synthetic_iid()
    adata2.obs["cat1"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    adata2.obs["cat2"] = np.random.randint(0, 5, size=(adata2.shape[0],))
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))

    artifact1_cov = ln.Artifact.from_anndata(adata1, key="part_one_cov.h5ad").save()
    artifact2_cov = ln.Artifact.from_anndata(adata2, key="part_two_cov.h5ad").save()

    collection_cov = ln.Collection([artifact1_cov, artifact2_cov], key="gather")
    collection_cov.save()

    artifacts_cov = collection_cov.artifacts.all()
    artifacts_cov.df()

    datamodule = MappedCollectionDataModule(
        collection_cov,
        batch_key="batch",
        batch_size=1024,
        join="inner",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    # pprint(datamodule.registry)

    model = scvi.model.SCVI(registry=datamodule.registry)
    # pprint(model.summary_stats)
    # pprint(model.module)

    model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    model.history.keys()

    # The way to extract the internal model analysis is by the inference_dataloader
    # Datamodule will always require passing it into all downstream functions.
    inference_dataloader = datamodule.inference_dataloader()
    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)
    _ = model.posterior_predictive_sample(dataloader=inference_dataloader)
    _ = model.get_normalized_expression(dataloader=inference_dataloader)
    _ = model.get_likelihood_parameters(dataloader=inference_dataloader)
    _ = model._get_denoised_samples(dataloader=inference_dataloader)
    _ = model.get_latent_library_size(dataloader=inference_dataloader, give_mean=False)

    # repeat but with other data with fewer indices and smaller batch size
    adata1_small = synthetic_iid(batch_size=10)
    adata1_small.obs["cat1"] = np.random.randint(0, 5, size=(adata1_small.shape[0],))
    adata1_small.obs["cat2"] = np.random.randint(0, 5, size=(adata1_small.shape[0],))
    adata1_small.obs["cont1"] = np.random.normal(size=(adata1_small.shape[0],))
    adata1_small.obs["cont2"] = np.random.normal(size=(adata1_small.shape[0],))
    adata2_small = synthetic_iid(batch_size=10)
    adata2_small.obs["cat1"] = np.random.randint(0, 5, size=(adata2_small.shape[0],))
    adata2_small.obs["cat2"] = np.random.randint(0, 5, size=(adata2_small.shape[0],))
    adata2_small.obs["cont1"] = np.random.normal(size=(adata2_small.shape[0],))
    adata2_small.obs["cont2"] = np.random.normal(size=(adata2_small.shape[0],))
    artifact1_small_cov = ln.Artifact.from_anndata(adata1_small, key="part_one_small.h5ad").save()
    artifact2_small_cov = ln.Artifact.from_anndata(adata2_small, key="part_two_small.h5ad").save()
    collection_small_cov = ln.Collection([artifact1_small_cov, artifact2_small_cov], key="gather")
    collection_small_cov.save()
    datamodule_small = MappedCollectionDataModule(
        collection_small_cov,
        batch_key="batch",
        batch_size=1024,
        join="inner",
        collection_val=collection_cov,
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
    )
    inference_dataloader_small = datamodule_small.inference_dataloader(batch_size=128)
    _ = model.get_elbo(return_mean=False, dataloader=inference_dataloader_small)
    _ = model.get_marginal_ll(n_mc_samples=3, dataloader=inference_dataloader_small)
    _ = model.get_reconstruction_error(return_mean=False, dataloader=inference_dataloader_small)
    _ = model.get_latent_representation(dataloader=inference_dataloader_small)
    _ = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["gene_1", "gene_2"], dataloader=inference_dataloader_small
    )
    _ = model.get_normalized_expression(n_samples=2, dataloader=inference_dataloader_small)

    # load and save and make a query with the other data
    model.save("lamin_model_cov", save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis
    loaded_model = scvi.model.SCVI.load("lamin_model_cov", adata=False)
    loaded_model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )


@pytest.mark.dataloader
@dependencies("tiledbsoma")
@dependencies("cellxgene_census")
def test_census_custom_dataloader_scvi(save_path: str):
    import cellxgene_census
    import tiledbsoma as soma

    # load census
    census = cellxgene_census.open_soma(census_version="2023-12-15")

    # do obs filtering (in this test we take a small dataset)
    experiment_name = "mus_musculus"
    obs_value_filter = (
        'is_primary_data == True and tissue_general in ["liver","heart"] and nnz >= 5000'
    )

    # to save time in this test, we manually filter var
    hv_idx = np.arange(10)  # just to make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use the class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_size=1024,
        shuffle=True,
        train_size=0.9,
        seed=42,
        batch_column_names=batch_keys,
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    n_layers = 1
    n_latent = 5

    # pprint(datamodule.registry)

    # creating the dataloader for trainset
    datamodule.setup()

    # We can now create the scVI model object and train it:
    model = scvi.model.SCVI(
        adata=None,
        registry=datamodule.registry,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood="nb",
        encode_covariates=False,
    )

    model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)
    model.history.keys()

    # save the model
    model.save("census_model", save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis
    scvi.model.SCVI.load("census_model", adata=False)

    # Generate cell embeddings
    inference_datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_labels=datamodule.batch_labels,
        batch_size=1024,
        shuffle=False,
        batch_column_names=batch_keys,
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    inference_datamodule.setup()

    # Datamodule will always require passing it into all downstream functions.
    # need to init the inference_dataloader before each of those commands:
    latent = model.get_latent_representation(
        dataloader=inference_datamodule.inference_dataloader()
    )
    print(latent.shape)
    _ = model.get_elbo(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_marginal_ll(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_reconstruction_error(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_latent_representation(dataloader=inference_datamodule.inference_dataloader())
    _ = model.posterior_predictive_sample(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_normalized_expression(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_likelihood_parameters(dataloader=inference_datamodule.inference_dataloader())
    _ = model._get_denoised_samples(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_latent_library_size(
        dataloader=inference_datamodule.inference_dataloader(), give_mean=False
    )

    # generating data from this census
    adata = cellxgene_census.get_anndata(
        census,
        organism=experiment_name,
        obs_value_filter=obs_value_filter,
        var_coords=hv_idx,
    )
    # verify cell order:
    assert np.array_equal(
        np.array(adata.obs["soma_joinid"]),
        inference_datamodule.train_dataset.query_ids.obs_joinids,
    )

    adata.obsm["scvi"] = latent

    # Additional things we would like to check
    # we make the batch name the same as in the model
    adata.obs["batch"] = adata.obs[batch_keys].agg("//".join, axis=1).astype("category")

    # query data
    model_query = model.load_query_data(
        adata=False, reference_model="census_model", registry=datamodule.registry
    )
    model_query.train(
        max_epochs=1, datamodule=datamodule, check_val_every_n_epoch=1, train_size=0.9
    )
    model_query.history.keys()

    _ = model_query.get_elbo(dataloader=inference_datamodule.inference_dataloader())
    _ = model_query.get_marginal_ll(dataloader=inference_datamodule.inference_dataloader())
    _ = model_query.get_reconstruction_error(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_query.get_latent_representation(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_query.posterior_predictive_sample(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_query.get_normalized_expression(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_query.get_likelihood_parameters(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_query._get_denoised_samples(dataloader=inference_datamodule.inference_dataloader())
    _ = model_query.get_latent_library_size(
        dataloader=inference_datamodule.inference_dataloader(), give_mean=False
    )

    scvi.model.SCVI.prepare_query_anndata(adata, "census_model", return_reference_var_names=True)
    scvi.model.SCVI.load_query_data(registry=datamodule.registry, reference_model="census_model")

    scvi.model.SCVI.prepare_query_anndata(adata, model)

    model.save("census_model2", save_anndata=False, overwrite=True, datamodule=datamodule)

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model_census3 = scvi.model.SCVI.load("census_model2", adata=adata)

    model_census3.train(
        datamodule=datamodule,
        max_epochs=1,
        check_val_every_n_epoch=1,
        train_size=0.9,
        early_stopping=False,
    )

    user_attributes_model_census3 = model_census3._get_user_attributes()
    pprint(user_attributes_model_census3)
    _ = model_census3.get_elbo()
    _ = model_census3.get_marginal_ll()
    _ = model_census3.get_reconstruction_error()
    _ = model_census3.get_latent_representation()
    _ = model_census3.posterior_predictive_sample()
    _ = model_census3.get_normalized_expression()
    _ = model_census3.get_likelihood_parameters()
    _ = model_census3._get_denoised_samples()
    _ = model_census3.get_latent_library_size(give_mean=False)
    for gene_likelihood in ["zinb", "nb", "poisson"]:
        model_adata = scvi.model.SCVI(adata, gene_likelihood=gene_likelihood)
        model_adata.train(1, check_val_every_n_epoch=1, train_size=0.9)

    scvi.model.SCVI.prepare_query_anndata(adata, "census_model2", return_reference_var_names=True)
    scvi.model.SCVI.load_query_data(adata, "census_model2")

    scvi.model.SCVI.prepare_query_anndata(adata, model_census3)
    scvi.model.SCVI.load_query_data(adata, model_census3)


@pytest.mark.dataloader
@dependencies("tiledbsoma")
@dependencies("cellxgene_census")
def test_census_custom_dataloader_scanvi(save_path: str):
    import cellxgene_census
    import tiledbsoma as soma

    # load census
    census = cellxgene_census.open_soma(census_version="2023-12-15")

    # do obs filtering (in this test we take a small dataset)
    experiment_name = "mus_musculus"
    obs_value_filter = (
        'is_primary_data == True and tissue_general in ["liver","heart"] and nnz >= 5000'
    )

    # to save time in this test, we manually filter var
    hv_idx = np.arange(10)  # just to make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use the class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    label_keys = ["tissue_general"]
    datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_size=1024,
        shuffle=True,
        seed=42,
        batch_column_names=batch_keys,
        label_keys=label_keys,
        train_size=0.9,
        model_name="SCANVI",
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    n_layers = 1
    n_latent = 5

    # pprint(datamodule.registry)

    # creating the dataloader for trainset
    datamodule.setup()

    # We can now create the scVI model object and train it:
    model = scvi.model.SCANVI(
        adata=None,
        registry=datamodule.registry,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood="nb",
        encode_covariates=False,
        datamodule=datamodule,
    )

    model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
        n_samples_per_label=100,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)
    model.history.keys()

    # save the model
    model.save("census_model_scanvi", save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis
    loaded_model = scvi.model.SCANVI.load(
        "census_model_scanvi", adata=False, datamodule=datamodule
    )
    loaded_model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
        n_samples_per_label=100,
    )

    # Generate cell embeddings
    inference_datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_labels=datamodule.batch_labels,
        batch_size=1024,
        shuffle=False,
        batch_column_names=batch_keys,
        label_keys=label_keys,
        model_name="SCANVI",
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    inference_datamodule.setup()

    model.predict(dataloader=inference_datamodule.inference_dataloader(), soft=False)

    latent = model.get_latent_representation(
        dataloader=inference_datamodule.inference_dataloader()
    )
    print(latent.shape)
    _ = model.get_elbo(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_marginal_ll(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_reconstruction_error(dataloader=inference_datamodule.inference_dataloader())
    _ = model.posterior_predictive_sample(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_normalized_expression(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_likelihood_parameters(dataloader=inference_datamodule.inference_dataloader())
    _ = model._get_denoised_samples(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_latent_library_size(
        dataloader=inference_datamodule.inference_dataloader(), give_mean=False
    )

    logged_keys = model.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "train_classification_loss" in logged_keys
    assert "train_accuracy" in logged_keys
    assert "train_f1_score" in logged_keys
    assert "train_calibration_error" in logged_keys
    assert "kl_global_validation" in logged_keys
    assert "kl_global_train" in logged_keys

    # train from scvi model
    model_scvi = scvi.model.SCVI(registry=datamodule.registry)

    model_scvi.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
        check_val_every_n_epoch=1,
        train_size=0.9,
    )
    model_scvi.save("census_model_scvi", save_anndata=False, overwrite=True, datamodule=datamodule)
    # We can now create the scVI model object and train it:
    model_scanvi_from_scvi = scvi.model.SCANVI.from_scvi_model(
        scvi_model=model_scvi,
        adata=None,
        registry=datamodule.registry,
        encode_covariates=False,
        unlabeled_category="label_0",
        labels_key="labels",
        datamodule=datamodule,
    )
    model_scanvi_from_scvi.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )
    # # save the model
    model_scanvi_from_scvi.save(
        "census_model_scanvi_from_scvi", save_anndata=False, overwrite=True, datamodule=datamodule
    )
    # # load it back and do downstream analysis
    model_scanvi_from_scvi_loaded = scvi.model.SCANVI.load(
        "census_model_scanvi_from_scvi", adata=False, datamodule=datamodule
    )
    model_scanvi_from_scvi_loaded.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )

    model_scanvi_from_scvi_loaded.predict(
        dataloader=inference_datamodule.inference_dataloader(), soft=False
    )

    latent = model_scanvi_from_scvi_loaded.get_latent_representation(
        dataloader=inference_datamodule.inference_dataloader()
    )
    print(latent.shape)
    _ = model_scanvi_from_scvi_loaded.get_elbo(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.get_marginal_ll(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.get_reconstruction_error(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.posterior_predictive_sample(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.get_normalized_expression(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.get_likelihood_parameters(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded._get_denoised_samples(
        dataloader=inference_datamodule.inference_dataloader()
    )
    _ = model_scanvi_from_scvi_loaded.get_latent_library_size(
        dataloader=inference_datamodule.inference_dataloader(), give_mean=False
    )

    logged_keys = model_scanvi_from_scvi_loaded.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "train_classification_loss" in logged_keys
    assert "train_accuracy" in logged_keys
    assert "train_f1_score" in logged_keys
    assert "train_calibration_error" in logged_keys
    assert "kl_global_validation" in logged_keys
    assert "kl_global_train" in logged_keys

    # generating adata from this census
    adata = cellxgene_census.get_anndata(
        census,
        organism=experiment_name,
        obs_value_filter=obs_value_filter,
        var_coords=hv_idx,
    )
    # verify cell order:
    assert np.array_equal(
        np.array(adata.obs["soma_joinid"]),
        inference_datamodule.train_dataset.query_ids.obs_joinids,
    )

    adata.obsm["scvi"] = latent

    # Additional things we would like to check
    # we make the batch name the same as in the model
    adata.obs["batch"] = adata.obs[batch_keys].agg("//".join, axis=1).astype("category")

    scvi.model.SCANVI.setup_anndata(
        adata, batch_key="batch", labels_key="tissue", unlabeled_category="label_0"
    )
    model_query_adata = model.load_query_data(
        adata=adata, reference_model="census_model_scanvi_from_scvi", datamodule=datamodule
    )
    model_query_adata.train(max_epochs=1, check_val_every_n_epoch=1, train_size=0.9)
    model_query_adata.predict(adata=adata)


@pytest.mark.dataloader
@dependencies("tiledbsoma")
@dependencies("cellxgene_census")
def test_census_custom_dataloader_scvi_with_covariates(save_path: str):
    import cellxgene_census
    import tiledbsoma as soma

    # load census
    census = cellxgene_census.open_soma(census_version="2023-12-15")

    # do obs filtering (in this test we take a small dataset)
    experiment_name = "mus_musculus"
    obs_value_filter = (
        'is_primary_data == True and tissue_general in ["liver","heart"] and nnz >= 5000'
    )

    # to save time in this test, we manually filter var
    hv_idx = np.arange(10)  # just to make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use the class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_size=1024,
        shuffle=True,
        train_size=0.9,
        seed=42,
        batch_column_names=batch_keys,
        categorical_covariate_keys=["sex"],
        continuous_covariate_keys=["raw_mean_nnz"],
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    n_layers = 1
    n_latent = 5

    # pprint(datamodule.registry)

    # creating the dataloader for trainset
    datamodule.setup()

    # We can now create the scVI model object and train it:
    model = scvi.model.SCVI(
        adata=None,
        registry=datamodule.registry,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood="nb",
        encode_covariates=False,
    )

    model.train(
        datamodule=datamodule,
        max_epochs=1,
        batch_size=1024,
        train_size=0.9,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)
    model.history.keys()

    # save the model
    model.save("census_model_cov", save_anndata=False, overwrite=True, datamodule=datamodule)

    # Generate cell embeddings
    inference_datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_labels=datamodule.batch_labels,
        batch_size=1024,
        shuffle=False,
        batch_column_names=batch_keys,
        categorical_covariate_keys=["sex"],
        continuous_covariate_keys=["raw_mean_nnz"],
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    inference_datamodule.setup()

    # Datamodule will always require passing it into all downstream functions.
    # need to init the inference_dataloader before each of those commands:
    latent = model.get_latent_representation(
        dataloader=inference_datamodule.inference_dataloader()
    )
    print(latent.shape)
    _ = model.get_elbo(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_marginal_ll(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_reconstruction_error(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_latent_representation(dataloader=inference_datamodule.inference_dataloader())
    _ = model.posterior_predictive_sample(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_normalized_expression(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_likelihood_parameters(dataloader=inference_datamodule.inference_dataloader())
    _ = model._get_denoised_samples(dataloader=inference_datamodule.inference_dataloader())
    _ = model.get_latent_library_size(
        dataloader=inference_datamodule.inference_dataloader(), give_mean=False
    )

    # generating data from this census
    adata = cellxgene_census.get_anndata(
        census,
        organism=experiment_name,
        obs_value_filter=obs_value_filter,
        var_coords=hv_idx,
    )
    # verify cell order:
    assert np.array_equal(
        np.array(adata.obs["soma_joinid"]),
        inference_datamodule.train_dataset.query_ids.obs_joinids,
    )

    adata.obsm["scvi"] = latent

    # Additional things we would like to check
    # we make the batch name the same as in the model
    adata.obs["batch"] = adata.obs[batch_keys].agg("//".join, axis=1).astype("category")

    model_census3 = scvi.model.SCVI.load("census_model_cov", adata=False)

    model_census3.train(
        datamodule=datamodule,
        max_epochs=1,
        early_stopping=False,
    )


@pytest.mark.dataloader
def test_annbatch(save_path: str):
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    path1 = os.path.join(save_path, "file1.h5ad")
    path2 = os.path.join(save_path, "file2.h5ad")
    collection_path = os.path.join(save_path, "annbatch_collection")

    adata1 = scvi.data.synthetic_iid(batch_size=20000)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=20000)
    adata2.X = csr_matrix(adata2.X)
    adata1.write(path1)
    adata2.write(path2)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        batch_size=4096,
        dataset_size=2_097_152,
    )

    assert dm.n_batch == 2, f"Expected 2 batches, got {dm.n_batch}"
    assert dm.n_labels > 0, f"Expected labels, got {dm.n_labels}"
    print(f"n_batch={dm.n_batch}, n_labels={dm.n_labels}, n_obs={dm.n_obs}, n_vars={dm.n_vars}")

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    model.history.keys()

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    print(latent.shape)

    # Validation: reuse the same zarr collection (collection_path_val == collection_path)
    dm_val = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        batch_size=4096,
        paths_val=[path1, path2],
        collection_path_val=collection_path,
    )

    model_val = scvi.model.SCVI(registry=dm_val.registry)
    model_val.train(
        max_epochs=1,
        datamodule=dm_val,
        check_val_every_n_epoch=1,
        train_size=0.9,
    )

    logged_keys = model_val.history.keys()
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "validation_loss" in logged_keys

    inference_dl_val = dm_val.inference_dataloader()
    _ = model_val.get_elbo(dataloader=inference_dl_val)
    _ = model_val.get_latent_representation(dataloader=inference_dl_val)


@pytest.mark.dataloader
def test_annbatch_with_covariates(save_path: str):
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    path1 = os.path.join(save_path, "file1.h5ad")
    path2 = os.path.join(save_path, "file2.h5ad")
    collection_path = os.path.join(save_path, "annbatch_covariates_collection")

    adata1 = scvi.data.synthetic_iid(batch_size=20000)
    adata1.X = csr_matrix(adata1.X)
    adata1.obs["cat1"] = np.random.randint(0, 5, size=(adata1.shape[0],)).astype(str)
    adata1.obs["cat2"] = np.random.randint(0, 3, size=(adata1.shape[0],)).astype(str)
    adata1.obs["cont1"] = np.random.normal(size=(adata1.shape[0],))
    adata1.obs["cont2"] = np.random.normal(size=(adata1.shape[0],))

    adata2 = scvi.data.synthetic_iid(batch_size=20000)
    adata2.X = csr_matrix(adata2.X)
    adata2.obs["cat1"] = np.random.randint(0, 5, size=(adata2.shape[0],)).astype(str)
    adata2.obs["cat2"] = np.random.randint(0, 3, size=(adata2.shape[0],)).astype(str)
    adata2.obs["cont1"] = np.random.normal(size=(adata2.shape[0],))
    adata2.obs["cont2"] = np.random.normal(size=(adata2.shape[0],))

    adata1.write(path1)
    adata2.write(path2)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
        batch_size=4096,
        dataset_size=2_097_152,
    )

    assert dm.n_batch == 2
    assert dm.n_labels > 0
    print(f"n_batch={dm.n_batch}, n_labels={dm.n_labels}")

    reg = dm.registry
    assert (
        reg["field_registries"]["extra_categorical_covs"]["summary_stats"][
            "n_extra_categorical_covs"
        ]
        == 2
    )
    assert (
        reg["field_registries"]["extra_continuous_covs"]["summary_stats"][
            "n_extra_continuous_covs"
        ]
        == 2
    )

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)

    logged_keys = model.history.keys()
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys

    inference_dl = dm.inference_dataloader()
    _ = model.get_elbo(dataloader=inference_dl)
    _ = model.get_latent_representation(dataloader=inference_dl)


@pytest.mark.dataloader
def test_annbatch_setup_scvi(save_path: str):
    """Test SCVI.setup_annbatch: build-once, reuse, and basic inference."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    collection_path = os.path.join(save_path, "annbatch_setup_scvi.zarr")

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)

    path1 = os.path.join(save_path, "setup_file1.h5ad")
    path2 = os.path.join(save_path, "setup_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    # --- First call: builds the zarr collection ---
    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        dataset_size=1024,
    )

    assert dm.n_batch == 2
    assert dm.n_vars == adata1.n_vars
    # Registry must carry actual gene names, not generic gene_i placeholders
    col_names = dm.registry["field_registries"]["X"]["state_registry"]["column_names"]
    assert col_names == list(adata1.var_names), "var_names not propagated into registry"

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history.keys()

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs
    assert latent.shape[1] == model.module.n_latent

    # --- Second call: reuses existing collection without re-specifying paths ---
    dm2 = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        batch_key="batch",
        batch_size=256,
    )
    assert dm2.n_batch == 2
    assert dm2.n_vars == adata1.n_vars

    # --- Third call: explicit rebuild (paths required) ---
    dm3 = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        rebuild=True,
    )
    model3 = scvi.model.SCVI(registry=dm3.registry)
    model3.train(max_epochs=1, datamodule=dm3)
    assert "elbo_train" in model3.history.keys()

    # --- With layer: save counts to a layer and load via layer= ---
    collection_layer_path = os.path.join(save_path, "annbatch_setup_scvi_layer.zarr")
    adata1_layer = adata1.copy()
    adata2_layer = adata2.copy()
    adata1_layer.layers["counts"] = adata1_layer.X.copy()
    adata2_layer.layers["counts"] = adata2_layer.X.copy()
    path1_layer = os.path.join(save_path, "setup_layer_file1.h5ad")
    path2_layer = os.path.join(save_path, "setup_layer_file2.h5ad")
    adata1_layer.write(path1_layer)
    adata2_layer.write(path2_layer)

    dm_layer = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_layer_path,
        paths=[path1_layer, path2_layer],
        batch_key="batch",
        layer="counts",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm_layer.n_vars == adata1.n_vars
    model_layer = scvi.model.SCVI(registry=dm_layer.registry)
    model_layer.train(max_epochs=1, datamodule=dm_layer)
    assert "elbo_train" in model_layer.history.keys()

    # --- With validation split ---
    collection_path_val = os.path.join(save_path, "annbatch_setup_scvi_val.zarr")
    dm_with_val = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        paths_val=[path1, path2],
        collection_path_val=collection_path_val,
    )
    model_val = scvi.model.SCVI(registry=dm_with_val.registry)
    model_val.train(max_epochs=1, datamodule=dm_with_val, check_val_every_n_epoch=1)


@pytest.mark.dataloader
def test_annbatch_setup_scanvi(save_path: str):
    """Test SCANVI.setup_annbatch: train, predict, and inference."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    collection_path = os.path.join(save_path, "annbatch_scanvi.zarr")

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)

    path1 = os.path.join(save_path, "scanvi_file1.h5ad")
    path2 = os.path.join(save_path, "scanvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    dm = scvi.model.SCANVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        unlabeled_category="Unknown",
        batch_size=256,
        dataset_size=1024,
    )

    assert dm.n_batch == 2
    assert dm.n_labels > 0
    assert dm.registry["model_name"] == "SCANVI"
    assert (
        dm.registry["field_registries"]["labels"]["state_registry"]["unlabeled_category"]
        == "Unknown"
    )

    model = scvi.model.SCANVI(
        adata=None,
        registry=dm.registry,
        encode_covariates=False,
        datamodule=dm,
    )
    model.train(max_epochs=1, datamodule=dm)

    logged_keys = model.history.keys()
    assert "elbo_train" in logged_keys
    assert "train_classification_loss" in logged_keys
    assert "train_accuracy" in logged_keys

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    predictions = model.predict(dataloader=inference_dl, soft=False)
    assert len(predictions) == dm.n_obs


@pytest.mark.dataloader
def test_annbatch_setup_base_sample_key(save_path: str):
    """Base setup_annbatch must forward sample_key to AnnbatchDataModule."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata1.obs["sample"] = "sample_A"
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    adata2.obs["sample"] = "sample_B"

    path1 = os.path.join(save_path, "base_sample_file1.h5ad")
    path2 = os.path.join(save_path, "base_sample_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "base_sample.zarr")
    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        sample_key="sample",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_samples == 2
    assert dm.registry["field_registries"]["sample"]["summary_stats"]["n_sample"] == 2


@pytest.mark.dataloader
def test_annbatch_setup_linear_scvi(save_path: str):
    """Test LinearSCVI.setup_annbatch: build, train, and inference."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "linear_scvi_file1.h5ad")
    path2 = os.path.join(save_path, "linear_scvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "linear_scvi.zarr")
    dm = scvi.model.LinearSCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_vars == adata1.n_vars

    model = scvi.model.LinearSCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs
    assert latent.shape[1] == model.n_latent


@pytest.mark.dataloader
def test_annbatch_setup_autozi(save_path: str):
    """Test AUTOZI.setup_annbatch: build, train."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "autozi_file1.h5ad")
    path2 = os.path.join(save_path, "autozi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "autozi.zarr")
    dm = scvi.model.AUTOZI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.model.AUTOZI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_condscvi(save_path: str):
    """Test CondSCVI.setup_annbatch: build, train."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "condscvi_file1.h5ad")
    path2 = os.path.join(save_path, "condscvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "condscvi.zarr")
    dm = scvi.model.CondSCVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_labels > 0

    model = scvi.model.CondSCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_decipher(save_path: str):
    """Test Decipher.setup_annbatch: build, train."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "decipher_file1.h5ad")
    path2 = os.path.join(save_path, "decipher_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "decipher.zarr")
    dm = scvi.external.Decipher.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_vars == adata1.n_vars

    model = scvi.external.Decipher(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_amortizedlda(save_path: str):
    """Test AmortizedLDA.setup_annbatch: build, train."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "lda_file1.h5ad")
    path2 = os.path.join(save_path, "lda_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "lda.zarr")
    dm = scvi.model.AmortizedLDA.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_vars == adata1.n_vars

    model = scvi.model.AmortizedLDA(registry=dm.registry, n_topics=5)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_sysvi(save_path: str):
    """Test SysVI.setup_annbatch: build, train with standard_normal prior."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "sysvi_file1.h5ad")
    path2 = os.path.join(save_path, "sysvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "sysvi.zarr")
    dm = scvi.external.SysVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.external.SysVI(registry=dm.registry, prior="standard_normal")
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_scar(save_path: str):
    """Test SCAR.setup_annbatch: build, train with provided ambient_profile."""
    import numpy as np
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    path1 = os.path.join(save_path, "scar_file1.h5ad")
    path2 = os.path.join(save_path, "scar_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "scar.zarr")
    dm = scvi.external.SCAR.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    # Provide ambient_profile as uniform prior (1/n_genes per gene)
    n_genes = adata1.n_vars
    ambient_profile = np.full((1, n_genes), 1.0 / n_genes, dtype=np.float32)

    model = scvi.external.SCAR(registry=dm.registry, ambient_profile=ambient_profile)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_mrvi(save_path: str):
    """Test MRVI/TorchMRVI.setup_annbatch with sample_key."""
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    adata1 = scvi.data.synthetic_iid(batch_size=500)
    adata1.X = csr_matrix(adata1.X)
    adata1.obs["donor"] = "donor_A"
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix(adata2.X)
    adata2.obs["donor"] = "donor_B"
    path1 = os.path.join(save_path, "mrvi_file1.h5ad")
    path2 = os.path.join(save_path, "mrvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "mrvi.zarr")
    dm = scvi.external.TorchMRVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        sample_key="donor",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_samples == 2
    assert dm.registry["field_registries"]["sample"]["summary_stats"]["n_sample"] == 2

    model = scvi.external.TorchMRVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_setup_peakvi(save_path: str):
    """Test PEAKVI.setup_annbatch: build, train with ATAC data."""
    import numpy as np
    import zarr
    from scipy.sparse import csr_matrix

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})

    # PEAKVI uses binary accessibility data
    adata1 = scvi.data.synthetic_iid(batch_size=500)
    # Convert to binary (0/1) to simulate ATAC
    adata1.X = csr_matrix((adata1.X > 0).astype(np.float32))
    adata2 = scvi.data.synthetic_iid(batch_size=500)
    adata2.X = csr_matrix((adata2.X > 0).astype(np.float32))
    path1 = os.path.join(save_path, "peakvi_file1.h5ad")
    path2 = os.path.join(save_path, "peakvi_file2.h5ad")
    adata1.write(path1)
    adata2.write(path2)

    collection_path = os.path.join(save_path, "peakvi.zarr")
    dm = scvi.model.PEAKVI.setup_annbatch(
        collection_path=collection_path,
        paths=[path1, path2],
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.model.PEAKVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm)
    assert "elbo_train" in model.history
