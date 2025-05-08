from __future__ import annotations

import os
from pprint import pprint

import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.dataloaders import MappedCollectionDataModule, TileDBDataModule
from scvi.utils import dependencies


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_small(save_path: str):
    os.system("lamin init --storage ./lamindb_collection")  # one time for github runner (comment)
    import lamindb as ln

    ln.setup.init()  # one time for github runner (comment out when runing localy)

    # prepare test data
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()

    artifact1 = ln.Artifact.from_anndata(adata1, key="part_one.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="part_two.h5ad").save()

    collection = ln.Collection([artifact1, artifact2], key="gather")
    collection.save()

    artifacts = collection.artifacts.all()
    artifacts.df()

    # large data example
    # ln.track("d1kl7wobCO1H0005")
    # ln.setup.init(name="lamindb_instance_name", storage=save_path)  # is this need in github test
    # ln.setup.init()
    # collection = ln.Collection.using("laminlabs/cellxgene").get(name="covid_normal_lung")
    # artifacts = collection.artifacts.all()
    # artifacts.df()

    datamodule = MappedCollectionDataModule(
        collection,
        batch_key="batch",
        batch_size=1024,
        join="inner",
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    pprint(datamodule.registry)

    model = scvi.model.SCVI(registry=datamodule.registry)
    pprint(model.summary_stats)
    pprint(model.module)

    model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    model.history.keys()

    # The way to extract the internal model analysis is by the inference_dataloader
    # Datamodule will always require to pass it into all downstream functions.
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

    # load and save and make query with the other data
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

    # cretae a regular model
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
def test_lamindb_dataloader_scanvi_small(save_path: str):
    # os.system("lamin init --storage ./lamindb_collection_scanvi") #(comment out runing localy)
    import lamindb as ln

    # ln.setup.init() # (comment out when runing localy)

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
def test_lamindb_dataloader_scvi_small_with_covariates(save_path: str):
    # os.system("lamin init --storage ./lamindb_collection_cov")  # one time for github runner
    import lamindb as ln

    # ln.setup.init()  # one time for github runner (comment out when runing localy)

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

    pprint(datamodule.registry)

    model = scvi.model.SCVI(registry=datamodule.registry)
    pprint(model.summary_stats)
    pprint(model.module)

    model.train(
        max_epochs=1,
        batch_size=1024,
        datamodule=datamodule,
    )
    model.history.keys()

    # The way to extract the internal model analysis is by the inference_dataloader
    # Datamodule will always require to pass it into all downstream functions.
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

    # load and save and make query with the other data
    model.save("lamin_model_cov", save_anndata=False, overwrite=True, datamodule=datamodule)


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

    # in order to save time in this test we manulay filter var
    hv_idx = np.arange(10)  # just to make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
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

    pprint(datamodule.registry)

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

    # Datamodule will always require to pass it into all downstream functions.
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

    # in order to save time in this test we manually filter var
    hv_idx = np.arange(10)  # just ot make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
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

    pprint(datamodule.registry)

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

    model.predict(dataloader=inference_datamodule.inference_dataloader(), soft=False)

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

    # in order to save time in this test we manulay filter var
    hv_idx = np.arange(10)  # just to make it smaller and faster for debug

    # For HVG, we can use the highly_variable_genes function provided in cellxgene_census,
    # which can compute HVGs in constant memory:
    hvg_query = census["census_data"][experiment_name].axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    )

    # We will now use class TileDBDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
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

    pprint(datamodule.registry)

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

    # Datamodule will always require to pass it into all downstream functions.
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
