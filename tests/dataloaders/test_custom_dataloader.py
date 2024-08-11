from __future__ import annotations

import cellxgene_census
import numpy as np
import pandas as pd
import pytest
import tiledbsoma as soma
from cellxgene_census.experimental.ml import experiment_dataloader
from cellxgene_census.experimental.ml.datamodule import CensusSCVIDataModule

import scvi
from scvi.data import synthetic_iid


@pytest.custom.dataloader
def test_custom_dataloader(save_path):
    # this test checks the local custom dataloder made by CZI and run several tests with it
    census = cellxgene_census.open_soma(census_version="stable")

    experiment_name = "mus_musculus"
    obs_value_filter = 'is_primary_data == True and tissue_general in ["kidney"] and nnz >= 3000'
    # top_n_hvg = 8000
    # hvg_batch = ["assay", "suspension_type"]
    #
    # # For HVG, we can use the `highly_variable_genes` function provided in the Census,
    # # which can compute HVGs in constant memory:
    #
    # query = census["census_data"][experiment_name].axis_query(
    #     measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=obs_value_filter)
    # )
    # hvgs_df = highly_variable_genes(query, n_top_genes=top_n_hvg, batch_key=hvg_batch)
    #
    # hv = hvgs_df.highly_variable
    # hv_idx = hv[hv].index

    hv_idx = np.arange(100)

    datamodule = CensusSCVIDataModule(
        census["census_data"][experiment_name],
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
        batch_size=1024,
        shuffle=True,
        batch_keys=["dataset_id", "assay", "suspension_type", "donor_id"],
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    datamodule.vars = hv_idx

    scvi.model._scvi.SCVI.setup_datamodule(datamodule)  # takes time

    adata = synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

    model = scvi.model.SCVI(adata, n_latent=10)
    model.train(max_epochs=1)

    dataloader = model._make_data_loader(adata)
    _ = model.get_elbo(dataloader=dataloader)
    _ = model.get_marginal_ll(dataloader=dataloader)
    _ = model.get_reconstruction_error(dataloader=dataloader)
    _ = model.get_latent_representation(dataloader=dataloader)

    scvi.model.SCVI.prepare_query_anndata(adata, reference_model=model)
    scvi.model.SCVI.load_query_data(adata, reference_model=model)

    n_layers = 1
    n_latent = 50

    model_census = scvi.model.SCVI(
        registry=datamodule.registry,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood="nb",
        encode_covariates=False,
    )

    batch_size = 1024
    train_size = 0.9
    max_epochs = 1

    model_census.train(
        datamodule=datamodule,
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=train_size,
        early_stopping=False,
    )

    user_attributes = model_census._get_user_attributes()
    user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

    model_census.save(save_path, overwrite=True)
    model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    model_census2.train(
        datamodule=datamodule,
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=train_size,
        early_stopping=False,
    )

    # takes time
    adata = cellxgene_census.get_anndata(
        census,
        organism=experiment_name,
        obs_value_filter=obs_value_filter,
        var_coords=hv_idx,
    )

    adata.obs["batch"] = (
        "batch_" + adata.obs[datamodule.batch_keys[0]].cat.codes.astype(str)
    ).astype("category")
    # adata.var_names = 'gene_'+adata.var_names #not sure we need it

    scvi.model.SCVI.prepare_query_anndata(adata, save_path)
    scvi.model.SCVI.load_query_data(registry=datamodule.registry, reference_model=save_path)

    scvi.model.SCVI.prepare_query_anndata(adata, model_census2)

    model_census3 = scvi.model.SCVI.load(save_path, adata=adata)

    scvi.model.SCVI.prepare_query_anndata(adata, save_path, return_reference_var_names=True)
    scvi.model.SCVI.load_query_data(adata, save_path)

    scvi.model.SCVI.prepare_query_anndata(adata, model_census3)
    scvi.model.SCVI.load_query_data(adata, model_census3)

    datamodule_inference = CensusSCVIDataModule(
        census["census_data"][experiment_name],
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
        batch_size=1024,
        shuffle=False,
        soma_chunk_size=50_000,
        batch_keys=["dataset_id", "assay", "suspension_type", "donor_id"],
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    datapipe = datamodule_inference.datapipe
    dataloader = experiment_dataloader(datapipe, num_workers=0, persistent_workers=False)
    mapped_dataloader = (
        datamodule_inference.on_before_batch_transfer(tensor, None) for tensor in dataloader
    )
    latent = model.get_latent_representation(dataloader=mapped_dataloader)

    emb_idx = datapipe._obs_joinids

    obs_soma_joinids = adata.obs["soma_joinid"]

    obs_indexer = pd.Index(emb_idx)
    idx = obs_indexer.get_indexer(obs_soma_joinids)
    # Reindexing is necessary to ensure that the cells in the embedding match the ones in
    # the anndata object.
    adata.obsm["scvi"] = latent[idx]
