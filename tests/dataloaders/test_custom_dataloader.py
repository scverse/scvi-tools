from __future__ import annotations

import os
from pprint import pprint
from time import time

import numpy as np
import torch

import scvi
from scvi.dataloaders import MappedCollectionDataModule, SCVIDataModule
from scvi.utils import dependencies


# @pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_scanvi(save_path: str = "."):
    os.system("lamin init --storage ./test-registries")
    import lamindb as ln

    # ln.setup.init(name="lamindb_instance_name", storage=save_path)  # is this need in github test
    # create synthetic data and than a collection of it

    collection = ln.Collection.get(name="3TNCsZZcnIBv2WGb0001")
    artifacts = collection.artifacts.all()
    artifacts.df()

    datamodule = MappedCollectionDataModule(
        collection, batch_key="assay", batch_size=1024, join="inner"
    )
    model = scvi.model.SCVI(adata=None, registry=datamodule.registry)
    pprint(model.summary_stats)
    pprint(model.module)
    inference_dataloader = datamodule.inference_dataloader()

    # Using regular adata laoder
    # adata = collection.load()  # try to compare this in regular settings
    # # setup large
    # SCVI.setup_anndata(
    #     adata,
    #     batch_key="assay",
    # )
    # model_reg = SCVI(adata)
    # start_time = time()
    # model_reg.train(max_epochs=10, batch_size=1024)
    # time_reg = time() - start_time
    # print(time_reg)

    start_time = time()
    model.train(max_epochs=10, batch_size=1024, datamodule=datamodule)
    time_lamin = time() - start_time
    print(time_lamin)

    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True)
    model_query = model.load_query_data(
        adata=False, reference_model="lamin_model", registry=datamodule.registry
    )
    model_query.train(max_epochs=1, datamodule=datamodule)
    _ = model_query.get_elbo(dataloader=inference_dataloader)
    _ = model_query.get_marginal_ll(dataloader=inference_dataloader)
    _ = model_query.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model_query.get_latent_representation(dataloader=inference_dataloader)

    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata=adata, reference_model="lamin_model")
    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata)
    model_query_adata.train(max_epochs=1)
    _ = model_query_adata.get_elbo()
    _ = model_query_adata.get_marginal_ll()
    _ = model_query_adata.get_reconstruction_error()
    _ = model_query_adata.get_latent_representation()
    _ = model_query_adata.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True)
    model.load("lamin_model", adata=False)
    model.load_query_data(adata=False, reference_model="lamin_model", registry=datamodule.registry)

    model.load_query_data(adata=adata, reference_model="lamin_model")
    model_adata = model.load("lamin_model", adata=adata)
    model_adata.train(max_epochs=1)
    model_adata.save("lamin_model_anndata", save_anndata=False, overwrite=True)
    model_adata.load("lamin_model_anndata", adata=False)
    model_adata.load_query_data(
        adata=False, reference_model="lamin_model_anndata", registry=datamodule.registry
    )


# @pytest.mark.dataloader
@dependencies("tiledbsoma")
@dependencies("cellxgene_census")
def test_census_custom_dataloader_scvi(save_path: str = "."):
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

    # We will now use class SCVIDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    datamodule = SCVIDataModule(
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

    # Setup the datamodule
    scvi.model._scvi.SCVI.setup_datamodule(datamodule)
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

    # what was on the notebook:
    # We can now save the trained model. As of the current writing, scvi-tools doesn't support
    # saving a model that wasn't generated through an AnnData loader, so we'll use some custom code
    model_state_dict = model.module.state_dict()
    var_names = hv_idx
    user_attributes = model._get_user_attributes()
    user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

    # save the model (TODO: Not working now in the regular way)
    model.save(save_path, save_anndata=False, overwrite=True)
    # load it back and do downstream analysis (not working)
    # model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    user_attributes.update(
        {
            "n_batch": datamodule.n_batch,
            "n_extra_categorical_covs": datamodule.registry["field_registries"][
                "extra_categorical_covs"
            ]["summary_stats"]["n_extra_categorical_covs"],
            "n_extra_continuous_covs": datamodule.registry["field_registries"][
                "extra_continuous_covs"
            ]["summary_stats"]["n_extra_continuous_covs"],
            "n_labels": datamodule.n_label,
            "n_vars": datamodule.n_vars,
            "batch_labels": datamodule.batch_labels,
        }
    )

    with open("model.pt", "wb") as f:
        torch.save(
            {
                "model_state_dict": model_state_dict,
                "var_names": var_names,
                "attr_dict": user_attributes,
            },
            f,
        )

    # We will now load the model back and use it to generate cell embeddings (the latent space)
    with open("model.pt", "rb") as f:
        torch_model = torch.load(f, weights_only=False)

        adict = torch_model["attr_dict"]
        params = adict["init_params_"]["non_kwargs"]

        n_batch = adict["n_batch"]
        # n_extra_categorical_covs = adict["n_extra_categorical_covs"]
        n_extra_continuous_covs = adict["n_extra_continuous_covs"]
        n_labels = adict["n_labels"]
        n_vars = adict["n_vars"]

        latent_distribution = params["latent_distribution"]
        dispersion = params["dispersion"]
        n_hidden = params["n_hidden"]
        dropout_rate = params["dropout_rate"]
        gene_likelihood = params["gene_likelihood"]

        model = scvi.model.SCVI(
            n_layers=params["n_layers"],
            n_latent=params["n_latent"],
            gene_likelihood=params["gene_likelihood"],
            encode_covariates=False,
        )

        module = model._module_cls(
            n_input=n_vars,
            n_batch=n_batch,
            n_labels=n_labels,
            n_continuous_cov=n_extra_continuous_covs,
            n_cats_per_cov=None,
            n_hidden=n_hidden,
            n_latent=n_latent,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
        )
        model.module = module

        model.module.load_state_dict(torch_model["model_state_dict"])

        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        model.to_device(device)
        model.module.eval()
        model.is_trained = True

    # Generate cell embeddings
    inference_datamodule = SCVIDataModule(
        hvg_query,
        layer_name="raw",
        batch_labels=adict["batch_labels"],
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
    # elbo = model.get_elbo(dataloader=inference_dataloader)
    # marginal_ll = model.get_marginal_ll(dataloader=inference_dataloader)
    # get_reconstruction_error = model.get_reconstruction_error(dataloader=inference_dataloader)

    # generating UMAP
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

    adata.obs["batch"] = adata.obs[batch_keys].agg("//".join, axis=1).astype("category")

    scvi.model.SCVI.prepare_query_anndata(adata, save_path, return_reference_var_names=True)
    # scvi.model.SCVI.load_query_data(registry=datamodule.registry, reference_model=save_path)

    # scvi.model.SCVI.prepare_query_anndata(adata, model)

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")  # needed?
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


# @pytest.mark.dataloader
@dependencies("tiledbsoma")
@dependencies("cellxgene_census")
def test_census_custom_dataloader_scanvi(save_path: str = "."):
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

    # We will now use class SCVIDataModule to connect TileDB-SOMA-ML with PyTorch Lightning.
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    label_keys = ["tissue_general"]
    datamodule = SCVIDataModule(
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

    # Setup the datamodule
    scvi.model.SCANVI.setup_datamodule(datamodule)
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

    # what was on the notebook:
    # We can now save the trained model. As of the current writing, scvi-tools doesn't support
    # saving a model that wasn't generated through an AnnData loader, so we'll use some custom code
    # model_state_dict = model.module.state_dict()
    # var_names = hv_idx
    user_attributes = model._get_user_attributes()
    user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

    # save the model (TODO: Not working now in the regular way)
    # model.save(save_path, save_anndata=False, overwrite=True)
    # load it back and do downstream analysis (not working)
    # model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    user_attributes.update(
        {
            "n_batch": datamodule.n_batch,
            "n_extra_categorical_covs": datamodule.registry["field_registries"][
                "extra_categorical_covs"
            ]["summary_stats"]["n_extra_categorical_covs"],
            "n_extra_continuous_covs": datamodule.registry["field_registries"][
                "extra_continuous_covs"
            ]["summary_stats"]["n_extra_continuous_covs"],
            "n_labels": datamodule.n_label,
            "n_vars": datamodule.n_vars,
            "batch_labels": datamodule.batch_labels,
        }
    )

    # Generate cell embeddings
    inference_datamodule = SCVIDataModule(
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
