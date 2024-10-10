from __future__ import annotations

import os
from pprint import pprint

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scvi
from scvi.data import synthetic_iid


@pytest.mark.custom_dataloader
def test_czi_custom_dataloader(save_path):
    # local branch with fix only for this test
    import sys

    # should be ready for importing the cloned branch on a remote machine that runs github action
    sys.path.insert(
        0,
        "/home/runner/work/scvi-tools/scvi-tools/"
        "cellxgene-census/api/python/cellxgene_census/src",
    )
    sys.path.insert(0, "src")
    import cellxgene_census
    import tiledbsoma as soma
    from cellxgene_census.experimental.ml import experiment_dataloader
    from cellxgene_census.experimental.ml.datamodule import CensusSCVIDataModule
    # from cellxgene_census.experimental.pp import highly_variable_genes

    # this test checks the local custom dataloder made by CZI and run several tests with it
    census = cellxgene_census.open_soma(census_version="stable")

    experiment_name = "mus_musculus"
    obs_value_filter = 'is_primary_data == True and tissue_general in ["kidney"] and nnz >= 3000'

    # This is under comments just to save time (selecting highly varkable genes):
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

    hv_idx = np.arange(100)  # just ot make it smaller and faster for debug

    # this is CZI part to be taken once all is ready
    batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
    datamodule = CensusSCVIDataModule(
        census["census_data"][experiment_name],
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
        batch_size=1024,
        shuffle=True,
        batch_keys=batch_keys,
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    # table of genes should be filtered by soma_joinid - but we should keep the encoded indexes
    # This is nice to have and might be uses in the downstream analysis
    # var_df = census["census_data"][experiment_name].ms["RNA"].var.read().concat().to_pandas()
    # var_df = var_df.loc[var_df.soma_joinid.isin(
    #     list(datamodule.datapipe.var_query.coords[0] if datamodule.datapipe.var_query is not None
    #          else range(datamodule.n_vars)))]

    # basicaly we should mimin everything below to any model census in scvi
    adata_orig = synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata_orig, batch_key="batch")

    model = scvi.model.SCVI(adata_orig, n_latent=10)
    model.train(max_epochs=1)

    # TODO: do we need to apply those functions to any census model as is?
    dataloader = model._make_data_loader(adata_orig)
    _ = model.get_elbo(dataloader=dataloader)
    _ = model.get_marginal_ll(dataloader=dataloader)
    _ = model.get_reconstruction_error(dataloader=dataloader)
    _ = model.get_latent_representation(dataloader=dataloader)

    scvi.model.SCVI.prepare_query_anndata(adata_orig, reference_model=model)
    scvi.model.SCVI.load_query_data(adata_orig, reference_model=model)

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)

    n_layers = 1
    n_latent = 50

    scvi.model._scvi.SCVI.setup_datamodule(datamodule)  # takes time
    model_census = scvi.model.SCVI(
        registry=datamodule.registry,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood="nb",
        encode_covariates=False,
    )

    pprint(datamodule.registry)

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

    user_attributes_model_census = model_census._get_user_attributes()
    # # TODO: do we need to put inside
    # user_attributes_model_census = \
    #     {a[0]: a[1] for a in user_attributes_model_census if a[0][-1] == "_"}
    pprint(user_attributes_model_census)
    # dataloader_census = model_census._make_data_loader(datamodule.datapipe)
    # # this casus errors
    # _ = model_census.get_elbo(dataloader=dataloader_census)
    # _ = model_census.get_marginal_ll(dataloader=dataloader_census)
    # _ = model_census.get_reconstruction_error(dataloader=dataloader_census)
    # _ = model_census.get_latent_representation(dataloader=dataloader_census)

    model_census.save(save_path, overwrite=True)
    model_census2 = scvi.model.SCVI.load(save_path, adata=False)

    model_census2.train(
        datamodule=datamodule,
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=train_size,
        early_stopping=False,
    )

    user_attributes_model_census2 = model_census2._get_user_attributes()
    pprint(user_attributes_model_census2)
    # dataloader_census2 = model_census2._make_data_loader()
    # this casus errors
    # _ = model_census2.get_elbo()
    # _ = model_census2.get_marginal_ll()
    # _ = model_census2.get_reconstruction_error()
    # _ = model_census2.get_latent_representation()

    # takes time
    adata = cellxgene_census.get_anndata(
        census,
        organism=experiment_name,
        obs_value_filter=obs_value_filter,
        var_coords=hv_idx,
    )

    # TODO: do we need to put inside (or is it alrady pre-made) - perhaps need to tell CZI
    adata.obs["batch"] = adata.obs[batch_keys].agg("".join, axis=1).astype("category")

    scvi.model.SCVI.prepare_query_anndata(adata, save_path)
    scvi.model.SCVI.load_query_data(registry=datamodule.registry, reference_model=save_path)

    scvi.model.SCVI.prepare_query_anndata(adata, model_census2)

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")  # needed?
    model_census3 = scvi.model.SCVI.load(save_path, adata=adata)

    model_census3.train(
        datamodule=datamodule,
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=train_size,
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

    # Create a dataloder of a CZI module
    datapipe = datamodule_inference.datapipe
    dataloader = experiment_dataloader(datapipe, num_workers=0, persistent_workers=False)
    mapped_dataloader = (
        datamodule_inference.on_before_batch_transfer(tensor, None) for tensor in dataloader
    )
    # _ = model_census.get_elbo(dataloader=mapped_dataloader)
    # _ = model_census.get_marginal_ll(dataloader=mapped_dataloader)
    # _ = model_census.get_reconstruction_error(dataloader=mapped_dataloader)
    latent = model_census.get_latent_representation(dataloader=mapped_dataloader)

    emb_idx = datapipe._obs_joinids

    obs_soma_joinids = adata.obs["soma_joinid"]

    obs_indexer = pd.Index(emb_idx)
    idx = obs_indexer.get_indexer(obs_soma_joinids)
    # Reindexing is necessary to ensure that the cells in the embedding match the ones in
    # the anndata object.
    adata.obsm["scvi"] = latent[idx]

    # #We can now generate the neighbors and the UMAP (tutorials)
    # sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi")
    # sc.tl.umap(adata, neighbors_key="scvi")
    # sc.pl.umap(adata, color="dataset_id", title="SCVI")
    #
    # sc.pl.umap(adata, color="tissue_general", title="SCVI")
    #
    # sc.pl.umap(adata, color="cell_type", title="SCVI")


@pytest.mark.custom_dataloader
def test_lamindb_custom_dataloader(save_path):
    # initialize a local lamin database
    os.system("lamin init --storage ~/scdataloader2 --schema bionty")
    # os.system("lamin close")
    # os.system("lamin load scdataloader")

    # local branch with fix only for this test
    import sys

    # should be ready for importing the cloned branch on a remote machine that runs github action
    sys.path.insert(
        0,
        "/home/runner/work/scvi-tools/scvi-tools/" "scDataLoader/",
    )
    sys.path.insert(0, "src")
    import lamindb as ln
    import tqdm
    from scdataloader import Collator, DataModule, SimpleAnnDataset

    # import bionty as bt
    # from scdataloader import utils
    from scdataloader.preprocess import (
        LaminPreprocessor,
        additional_postprocess,
        additional_preprocess,
    )

    # import numpy as np
    # import tiledbsoma as soma
    from scdataloader.utils import populate_my_ontology
    from torch.utils.data import DataLoader
    # from scdataloader.base import NAME
    # from cellxgene_census.experimental.ml import experiment_dataloader

    # populate_my_ontology() #to populate everything (recommended) (can take 2-10mns)

    populate_my_ontology(
        organisms=["NCBITaxon:10090", "NCBITaxon:9606"],
        sex=["PATO:0000384", "PATO:0000383"],
    )

    # preprocess datasets - do we need this part?
    DESCRIPTION = "preprocessed by scDataLoader"

    cx_dataset = (
        ln.Collection.using(instance="laminlabs/cellxgene")
        .filter(name="cellxgene-census", version="2023-12-15")
        .one()
    )
    cx_dataset, len(cx_dataset.artifacts.all())

    do_preprocess = LaminPreprocessor(
        additional_postprocess=additional_postprocess,
        additional_preprocess=additional_preprocess,
        skip_validate=True,
        subset_hvg=0,
    )

    do_preprocess(cx_dataset, name=DESCRIPTION, description=DESCRIPTION, start_at=1, version="2")

    # create dataloaders

    datamodule = DataModule(
        collection_name="preprocessed dataset",
        organisms=["NCBITaxon:9606"],  # organism that we will work on
        how="most expr",  # for the collator (most expr genes only will be selected) / "some"
        max_len=1000,  # only the 1000 most expressed
        batch_size=64,
        num_workers=1,
        validation_split=0.1,
        test_split=0,
    )

    # we setup the datamodule (as exemplified in lightning's good practices, b
    # ut there might be some things to improve here)
    # testfiles = datamodule.setup()

    for i in tqdm.tqdm(datamodule.train_dataloader()):
        # pass #or do pass
        print(i)
        break

    # with lightning:
    # Trainer(model, datamodule)

    # Read adata and create lamindb dataloader
    adata_orig = sc.read_h5ad("./scDataLoader/tests/test.h5ad")
    # preprocessor = Preprocessor(do_postp=False)
    # adata = preprocessor(adata_orig)
    adataset = SimpleAnnDataset(adata_orig, obs_to_output=["organism_ontology_term_id"])
    col = Collator(
        organisms=["NCBITaxon:9606"],
        max_len=1000,
        how="random expr",
    )
    dataloader = DataLoader(
        adataset,
        collate_fn=col,
        batch_size=4,
        num_workers=1,
        shuffle=False,
    )

    # We will now create the SCVI model object:
    # Its parameters:
    # n_layers = 1
    # n_latent = 10
    # batch_size = 1024
    # train_size = 0.9
    # max_epochs = 1

    # def on_before_batch_transfer(
    #     batch: tuple[torch.Tensor, torch.Tensor],
    # ) -> dict[str, torch.Tensor | None]:
    #     """Format the datapipe output with registry keys for scvi-tools."""
    #     X, obs = batch
    #     X_KEY: str = "X"
    #     BATCH_KEY: str = "batch"
    #     LABELS_KEY: str = "labels"
    #     return {
    #         X_KEY: X,
    #         BATCH_KEY: obs,
    #         LABELS_KEY: None,
    #     }

    # Try the lamindb dataloder on a trained scvi-model with adata
    # adata = adata.copy()
    scvi.model.SCVI.setup_anndata(adata_orig, batch_key="cell_type_ontology_term_id")
    model = scvi.model.SCVI(adata_orig, n_latent=10)
    model.train(max_epochs=1)
    # dataloader2 = experiment_dataloader(dataloader, num_workers=0, persistent_workers=False)
    # mapped_dataloader = (
    #     on_before_batch_transfer(tensor, None) for tensor in dataloader2
    # )
    # dataloader = model._make_data_loader(mapped_dataloader)
    _ = model.get_elbo(dataloader=dataloader)
    _ = model.get_marginal_ll(dataloader=dataloader)
    _ = model.get_reconstruction_error(dataloader=dataloader)
    _ = model.get_latent_representation(dataloader=dataloader)

    # scvi.model._scvi.SCVI.setup_datamodule(datamodule)  # takes time
    # model_lamindb = scvi.model.SCVI(
    #     registry=datamodule.registry,
    #     n_layers=n_layers,
    #     n_latent=n_latent,
    #     gene_likelihood="nb",
    #     encode_covariates=False,
    # )
    #
    # pprint(datamodule.registry)
    #
    # model_lamindb.train(
    #     datamodule=datamodule,
    #     max_epochs=max_epochs,
    #     batch_size=batch_size,
    #     train_size=train_size,
    #     early_stopping=False,
    # )
    # We have to create a registry without setup_anndata that contains the same elements
    # The other way will be to fill the model ,LIKE IN CELLXGENE NOTEBOOK
    # need to pass here new object of registry taht contains everything we will need
