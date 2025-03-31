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
    os.system("lamin init --storage ./lamindb_collection")  # one time for github runner
    import lamindb as ln

    # from scipy.sparse import csc_matrix, csr_matrix
    # import dask
    # import spatialdata
    ln.setup.init()  # one time for github runner

    # prepare test data
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()

    artifact1 = ln.Artifact.from_anndata(adata1, key="part_one.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="part_two.h5ad").save()

    collection = ln.Collection([artifact1, artifact2], key="gather")
    # test mapped without saving first
    # with collection.mapped() as ls_ds:
    #    assert ls_ds.__class__.__name__ == "MappedCollection"
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
        collection, batch_key="batch", batch_size=1024, join="inner"
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    pprint(datamodule.registry)

    model = scvi.model.SCVI(adata=None, registry=datamodule.registry)
    pprint(model.summary_stats)
    pprint(model.module)
    inference_dataloader = datamodule.inference_dataloader()

    model.train(max_epochs=1, batch_size=1024, datamodule=datamodule)

    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True, datamodule=datamodule)
    model_query = model.load_query_data(
        adata=False, reference_model="lamin_model", registry=datamodule.registry
    )
    model_query.train(max_epochs=1, datamodule=datamodule)
    _ = model_query.get_elbo(dataloader=inference_dataloader)
    _ = model_query.get_marginal_ll(dataloader=inference_dataloader)
    _ = model_query.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model_query.get_latent_representation(dataloader=inference_dataloader)

    adata = collection.load(join="inner")
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    with pytest.raises(ValueError):
        model.load_query_data(adata=adata)
    model_query_adata = model.load_query_data(adata=adata, reference_model="lamin_model")
    model_query_adata.train(max_epochs=1)
    _ = model_query_adata.get_elbo()
    _ = model_query_adata.get_marginal_ll()
    _ = model_query_adata.get_reconstruction_error()
    _ = model_query_adata.get_latent_representation()
    _ = model_query_adata.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True, datamodule=datamodule)
    model.load("lamin_model", adata=False)
    model.load_query_data(adata=False, reference_model="lamin_model", registry=datamodule.registry)

    model.load_query_data(adata=adata, reference_model="lamin_model")
    model_adata = model.load("lamin_model", adata=adata)
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model_adata.train(max_epochs=1)
    model_adata.save(
        "lamin_model_anndata", save_anndata=True, overwrite=True, datamodule=datamodule
    )
    model_adata.load("lamin_model_anndata")
    model_adata.load_query_data(
        adata=adata, reference_model="lamin_model_anndata", registry=datamodule.registry
    )


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scanvi_small(save_path: str):
    # os.system("lamin init --storage ./lamindb_collection")
    import lamindb as ln

    # from scipy.sparse import csc_matrix, csr_matrix
    # import dask
    # import spatialdata
    # ln.setup.init()

    # prepare test data
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()

    artifact1 = ln.Artifact.from_anndata(adata1, key="part_one.h5ad").save()
    artifact2 = ln.Artifact.from_anndata(adata2, key="part_two.h5ad").save()

    collection = ln.Collection([artifact1, artifact2], key="gather")
    # test mapped without saving first
    # with collection.mapped() as ls_ds:
    #    assert ls_ds.__class__.__name__ == "MappedCollection"
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
        label_key="labels",
        batch_key="batch",
        batch_size=1024,
        join="inner",
        unlabeled_category="label_0",
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
        check_val_every_n_epoch=1,
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)

    # save the model
    # model.save(save_path, save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis (not working)
    # model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    inference_dataloader = datamodule.inference_dataloader()

    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)

    logged_keys = model.history.keys()
    # assert "elbo_validation" in logged_keys
    # assert "reconstruction_loss_validation" in logged_keys
    # assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    # assert "validation_classification_loss" in logged_keys
    # assert "validation_accuracy" in logged_keys
    # assert "validation_f1_score" in logged_keys
    # assert "validation_calibration_error" in logged_keys


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
        seed=42,
        batch_column_names=batch_keys,
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    n_layers = 1
    n_latent = 5

    pprint(datamodule.registry)

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
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)

    # save the model
    model.save(save_path, save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis (not working)
    scvi.model.SCVI.load(save_path, adata=False)

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

    # creating the dataloader
    inference_dataloader = (
        inference_datamodule.on_before_batch_transfer(batch, None)
        for batch in inference_datamodule.train_dataloader()
    )

    latent = model.get_latent_representation(dataloader=inference_dataloader)
    print(latent.shape)
    # need to init the inference_dataloader before each of those commands:
    # elbo = model.get_elbo(dataloader=inference_dataloader)
    # marginal_ll = model.get_marginal_ll(dataloader=inference_dataloader)
    # get_reconstruction_error = model.get_reconstruction_error(dataloader=inference_dataloader)

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
    # we cmake the batch name the same as in the model
    adata.obs["batch"] = adata.obs[batch_keys].agg("//".join, axis=1).astype("category")

    scvi.model.SCVI.prepare_query_anndata(adata, save_path, return_reference_var_names=True)
    # scvi.model.SCVI.load_query_data(registry=datamodule.registry, reference_model=save_path)

    scvi.model.SCVI.prepare_query_anndata(adata, model)

    model.save(save_path, save_anndata=False, overwrite=True, datamodule=datamodule)

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model_census3 = scvi.model.SCVI.load(save_path, adata=adata)

    model_census3.train(
        datamodule=datamodule,
        max_epochs=1,
        early_stopping=False,
    )

    user_attributes_model_census3 = model_census3._get_user_attributes()
    pprint(user_attributes_model_census3)
    _ = model_census3.get_elbo()
    _ = model_census3.get_marginal_ll()
    _ = model_census3.get_reconstruction_error()
    _ = model_census3.get_latent_representation()

    scvi.model.SCVI.prepare_query_anndata(adata, save_path, return_reference_var_names=True)
    scvi.model.SCVI.load_query_data(adata, save_path)

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

    # in order to save time in this test we manulay filter var
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
        unlabeled_category="label_0",
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    print(datamodule.n_obs, datamodule.n_vars, datamodule.n_batch)

    n_layers = 1
    n_latent = 5

    pprint(datamodule.registry)

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
        check_val_every_n_epoch=1,
        early_stopping=False,
    )

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)

    # save the model
    # model.save(save_path, save_anndata=False, overwrite=True, datamodule=datamodule)
    # load it back and do downstream analysis (not working)
    # model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    # Generate cell embeddings
    inference_datamodule = TileDBDataModule(
        hvg_query,
        layer_name="raw",
        batch_labels=datamodule.batch_labels,
        batch_size=1024,
        shuffle=False,
        batch_column_names=batch_keys,
        label_keys=label_keys,
        unlabeled_category="label_0",
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    inference_datamodule.setup()

    # creating the dataloader
    inference_dataloader = (
        inference_datamodule.on_before_batch_transfer(batch, None)
        for batch in inference_datamodule.train_dataloader()
    )

    latent = model.get_latent_representation(dataloader=inference_dataloader)
    print(latent.shape)
    # elbo = model.get_elbo(dataloader=inference_dataloader)
    # marginal_ll = model.get_marginal_ll(dataloader=inference_dataloader)
    # get_reconstruction_error = model.get_reconstruction_error(dataloader=inference_dataloader)

    logged_keys = model.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    # assert "validation_classification_loss" in logged_keys
    # assert "validation_accuracy" in logged_keys
    # assert "validation_f1_score" in logged_keys
    # assert "validation_calibration_error" in logged_keys
