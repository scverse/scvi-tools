# from __future__ import annotations
#
# import sys
#
# sys.path.insert(0, "/Users/orikr/Documents/cellxgene-census/api/python/cellxgene_census/src")
# sys.path.insert(0, "src")
#
# import cellxgene_census
# import numpy as np
# import tiledbsoma as soma
# from cellxgene_census.experimental.ml.datamodule import (
#     CensusSCVIDataModule,  # WE RAN FROM LOCAL LIB
# )
# from cellxgene_census.experimental.pp import highly_variable_genes
#
# import scvi
# from scvi.data import synthetic_iid
# from scvi.model import SCVI
#
# # cellxgene_census.__file__, scvi.__file__
#
# # We will now create the SCVI model object:
# # Its parameters:
# n_layers = 1
# n_latent = 10
# batch_size = 1024
# train_size = 0.9
# max_epochs = 1
#
# # We have to create a registry without setup_anndata that contains the same elements
# # The other way will be to fill the model ,LIKE IN CELLXGENE NOTEBOOK
# # need to pass here new object of registry taht contains everything we will need
#
# # First lets see CELLXGENE example using pytorch loaders implemented now in our repo
# census = cellxgene_census.open_soma(census_version="stable")
# experiment_name = "mus_musculus"
# obs_value_filter = 'is_primary_data == True and tissue_general in ["spleen"] and nnz >= 300'
# top_n_hvg = 800
# hvg_batch = ["assay", "suspension_type"]
#
# # THIS WILL TAKE FEW MINUTES TO RUN!
# query = census["census_data"][experiment_name].axis_query(
#     measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=obs_value_filter)
# )
# hvgs_df = highly_variable_genes(query, n_top_genes=top_n_hvg, batch_key=hvg_batch)
# hv = hvgs_df.highly_variable
# hv_idx = hv[hv].index
# hv_idx = np.arange(100)  # just randomly select smaller number of indices
#
# # Now load the custom data module CZI did that now exists in our db
# # (and we will later want to elaborate with more info from our original anndata registry)
# # This thing is done by the user in any form they want
# datamodule = CensusSCVIDataModule(
#     census["census_data"][experiment_name],
#     measurement_name="RNA",
#     X_name="raw",
#     obs_query=soma.AxisQuery(value_filter=obs_value_filter),
#     var_query=soma.AxisQuery(coords=(list(hv_idx),)),
#     batch_size=1024,
#     shuffle=True,
#     batch_keys=["dataset_id", "assay", "suspension_type", "donor_id"],
#     dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
# )
#
# datamodule.vars = hv_idx
#
#
# # The next part is the same as test_scvi_train_custom_dataloader
# def test_scvi_train_custom_dataloader(n_latent: int = 10):
#     adata = synthetic_iid()
#     scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
#     model = scvi.model.SCVI(adata, n_latent=n_latent)
#     model.train(max_epochs=1)
#     dataloader = model._make_data_loader(adata)
#     _ = model.get_elbo(dataloader=dataloader)
#     _ = model.get_marginal_ll(dataloader=dataloader)
#     _ = model.get_reconstruction_error(dataloader=dataloader)
#     _ = model.get_latent_representation(dataloader=dataloader)
#
#     scvi.model.SCVI.prepare_query_anndata(adata, model)
#     query_model = scvi.model.SCVI.load_query_data(adata, model)
#
#
# def test_scvi_train_custom_datamodule(datamodule=datamodule):
#     # This is a new func to implement
#     # will take a bit of time to end
#     SCVI.setup_datamodule(datamodule)
#
#     # We will now create the SCVI model object from custom data module:
#     model_census = scvi.model.SCVI(
#         registry=datamodule.registry,
#         n_layers=n_layers,
#         n_latent=n_latent,
#         gene_likelihood="nb",
#         encode_covariates=False,
#     )
#
#     # The CZI data module is a refined data module while SCVI is a lighting datamodule
#     # Altough this is only 1 epoch it will take few mins on local machine
#     model_census.train(
#         datamodule=datamodule,
#         max_epochs=max_epochs,
#         batch_size=batch_size,
#         train_size=train_size,
#         early_stopping=False,
#     )
#
#     # We can now save the trained model. As of the current writing date (June 2024),
#     # scvi-tools doesn't support saving a model that wasn't generated through an AnnData loader,
#     # so we'll use some custom code:
#     # model_state_dict = model_census.module.state_dict()
#     # var_names = hv_idx.to_numpy()
#     # user_attributes = model_census._get_user_attributes()
#     # user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}
#     model_census.save("dataloader_model2", overwrite=True)
#     model_census_loaded = scvi.model.SCVI.load("dataloader_model2", adata=False)
#
#
# def test_scvi_train_custom_datamodule_from_loaded_model(datamodule=datamodule):
#     model_census_loaded = scvi.model.SCVI.load("dataloader_model2", adata=False)
#
#     # see if can train from loaded models
#     model_census_loaded.train(
#         datamodule=datamodule,
#         max_epochs=max_epochs,
#         batch_size=batch_size,
#         train_size=train_size,
#         early_stopping=False,
#     )
#
#
# def test_scvi_get_anndata_load_anndatamodule_from_custom_datamodule(datamodule=datamodule):
#     # we will perform here several task that deals with transforming custom data module to
#     # the regular ann datamodule - will take a bit of time
#     adata = cellxgene_census.get_anndata(
#         census, organism=experiment_name, obs_value_filter=obs_value_filter, var_coords=hv_idx
#     )
#
#     adata = adata[:, datamodule.vars].copy()
#
#     adata.obs.head()
#
#     # ORI Replace this with the function to generate batch key used in the datamodule.
#     # "12967895-3d58-4e93-be2c-4e1bcf4388d510x 5' v1cellHCA_Mou_3"
#     adata.obs["batch"] = (
#         "batch_" + adata.obs[datamodule.batch_keys[0]].cat.codes.astype(str)
#     ).astype("category")
#     # adata.var_names = 'gene_'+adata.var_names #not sure we need it
#
#     # We will now load the model back and use it to generate cell embeddings (the latent space),
#     # which can then be used for further analysis. Note that we still need to use some custom
#     # code for loading the model, which includes loading the parameters from the `attr_dict` node
#     # stored in the model.
#
#     # loading and setupanndata
#     model_census2 = scvi.model.SCVI.load("dataloader_model2", adata=False)
#     model_census2.setup_anndata(adata, batch_key="batch")
#     # model_census2.adata = deepcopy(adata)
#
#     # ORI Works when loading from disk
#     scvi.model.SCVI.prepare_query_anndata(
#         adata, "dataloader_model2", return_reference_var_names=True
#     )
#     # ORI This one still needs to be fixed.
#     scvi.model.SCVI.prepare_query_anndata(adata, model_census2, return_reference_var_names=True)
#     query_model = scvi.model.SCVI.load_query_data(
#         registry=datamodule.registry, reference_model="dataloader_model2"
#     )
#
#     # ORI Should work when setting up the AnnData correctly. scANVI with DataModule is not yet
#     # supported as DataModule can't take a labels_key.
#     scanvae = scvi.model.SCANVI.from_scvi_model(
#         model_census2,
#         adata=adata,
#         unlabeled_category="Unknown",
#         labels_key="cell_type",
#     )
#
#     # ORI - check it should work with a model initialized with AnnData.
#     # See below not fully working yet
#     model_census3 = scvi.model.SCVI.load("dataloader_model2", adata=adata)
#
#     scvi.model.SCVI.prepare_query_anndata(
#         adata, "dataloader_model2", return_reference_var_names=True
#     )
#     query_model = scvi.model.SCVI.load_query_data(adata, "dataloader_model2")
#
#     scvi.model.SCVI.prepare_query_anndata(adata, model_census3)
#     query_model = scvi.model.SCVI.load_query_data(adata, model_census3)
