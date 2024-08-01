from __future__ import annotations

import sys

sys.path.insert(0, "/Users/orikr/Documents/cellxgene-census/api/python/cellxgene_census/src")
sys.path.insert(0, "src")

import cellxgene_census
import numpy as np
import tiledbsoma as soma
from cellxgene_census.experimental.ml.datamodule import (
    CensusSCVIDataModule,  # WE RAN FROM LOCAL LIB
)
from cellxgene_census.experimental.pp import highly_variable_genes

import scvi
from scvi.data import _constants, synthetic_iid
from scvi.utils import attrdict

# cellxgene_census.__file__, scvi.__file__

# We will now create the SCVI model object:
# Its parameters:
n_layers = 1
n_latent = 10
batch_size = 1024
train_size = 0.9
max_epochs = 1

# We have to create a registry without setup_anndata that contains the same elements
# The other way will be to fill the model ,LIKE IN CELLXGENE NOTEBOOK
# need to pass here new object of registry taht contains everything we will need

# First lets see CELLXGENE example using pytorch loaders implemented now in our repo
census = cellxgene_census.open_soma(census_version="stable")
experiment_name = "mus_musculus"
obs_value_filter = 'is_primary_data == True and tissue_general in ["spleen"] and nnz >= 300'
top_n_hvg = 8000
hvg_batch = ["assay", "suspension_type"]
# THIS WILL TAKE FEW MINUTES TO RUN!
query = census["census_data"][experiment_name].axis_query(
    measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=obs_value_filter)
)
hvgs_df = highly_variable_genes(query, n_top_genes=top_n_hvg, batch_key=hvg_batch)
hv = hvgs_df.highly_variable
hv_idx = hv[hv].index

# Now load the custom data module CZI did that now exists in our db
# (and we will later want to elaborate with more info from our original anndata registry)
# This thing is done by the user in any form they want
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


def _get_summary_stats_from_registry(registry: dict) -> attrdict:
    summary_stats = {}
    for field_registry in registry[_constants._FIELD_REGISTRIES_KEY].values():
        field_summary_stats = field_registry[_constants._SUMMARY_STATS_KEY]
        summary_stats.update(field_summary_stats)
    return attrdict(summary_stats)


def setup_datamodule(datamodule: CensusSCVIDataModule):
    datamodule.registry = {
        "scvi_version": scvi.__version__,
        "model_name": "SCVI",
        "setup_args": {
            "layer": None,
            "batch_key": "batch",
            "labels_key": None,
            "size_factor_key": None,
            "categorical_covariate_keys": None,
            "continuous_covariate_keys": None,
        },
        "field_registries": {
            "X": {
                "data_registry": {"attr_name": "X", "attr_key": None},
                "state_registry": {
                    "n_obs": datamodule.n_obs,
                    "n_vars": datamodule.n_vars,
                    "column_names": datamodule.vars,
                },
                "summary_stats": {"n_vars": datamodule.n_vars, "n_cells": datamodule.n_obs},
            },
            "batch": {
                "data_registry": {"attr_name": "obs", "attr_key": "_scvi_batch"},
                "state_registry": {
                    "categorical_mapping": datamodule.datapipe.obs_encoders["batch"].classes_,
                    "original_key": "batch",
                },
                "summary_stats": {"n_batch": datamodule.n_batch},
            },
            "labels": {
                "data_registry": {"attr_name": "obs", "attr_key": "_scvi_labels"},
                "state_registry": {
                    "categorical_mapping": np.array([0]),
                    "original_key": "_scvi_labels",
                },
                "summary_stats": {"n_labels": 1},
            },
            "size_factor": {"data_registry": {}, "state_registry": {}, "summary_stats": {}},
            "extra_categorical_covs": {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_categorical_covs": 0},
            },
            "extra_continuous_covs": {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_continuous_covs": 0},
            },
        },
        "setup_method_name": "setup_datamodule",
    }
    datamodule.summary_stats = _get_summary_stats_from_registry(datamodule.registry)
    datamodule.var_names = [str(i) for i in datamodule.vars]


# This is a new func to implement (Implemented Above but we need in our code base as well)
# will take a bit of time to end
setup_datamodule(datamodule)

# The next part is the same as test_scvi_train_custom_dataloader

adata = synthetic_iid()
scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=10)
model.train(max_epochs=1)
dataloader = model._make_data_loader(adata)
_ = model.get_elbo(dataloader=dataloader)
_ = model.get_marginal_ll(dataloader=dataloader)
_ = model.get_reconstruction_error(dataloader=dataloader)
_ = model.get_latent_representation(dataloader=dataloader)

# ORI I broke the code here also for standard models. Please first fix this. - it is fixed
scvi.model.SCVI.prepare_query_anndata(adata, model)
query_model = scvi.model.SCVI.load_query_data(adata, model)

# We will now create the SCVI model object:
model_census = scvi.model.SCVI(
    datamodule=datamodule,
    n_layers=n_layers,
    n_latent=n_latent,
    gene_likelihood="nb",
    encode_covariates=False,
)

# The CZI data module is a refined data module while SCVI is a lighting datamodule
# Altough this is only 1 epoch it will take few mins on local machine
model_census.train(
    datamodule=datamodule,
    max_epochs=max_epochs,
    batch_size=batch_size,
    train_size=train_size,
    early_stopping=False,
)

# We can now save the trained model. As of the current writing date (June 2024),
# scvi-tools doesn't support saving a model that wasn't generated through an AnnData loader,
# so we'll use some custom code:
# model_state_dict = model_census.module.state_dict()
# var_names = hv_idx.to_numpy()
# user_attributes = model_census._get_user_attributes()
# user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}
model_census.save("dataloader_model2", overwrite=True)

# We are now turning this data module back to AnnData
adata = cellxgene_census.get_anndata(
    census,
    organism=experiment_name,
    obs_value_filter=obs_value_filter,
)

adata = adata[:, datamodule.vars].copy()

adata.obs.head()

# ORI Replace this with the function to generate batch key used in the datamodule.
# "12967895-3d58-4e93-be2c-4e1bcf4388d510x 5' v1cellHCA_Mou_3"
adata.obs["batch"] = ("batch_" + adata.obs[datamodule.batch_keys[0]].cat.codes.astype(str)).astype(
    "category"
)
# adata.var_names = 'gene_'+adata.var_names #not sure we need it

# We will now load the model back and use it to generate cell embeddings (the latent space),
# which can then be used for further analysis. Note that we still need to use some custom code for
# loading the model, which includes loading the parameters from the `attr_dict` node stored in
# the model.

model_census2 = scvi.model.SCVI.load("dataloader_model2", datamodule=datamodule)
model_census2.setup_anndata(adata, batch_key="batch")
# model_census2.adata = deepcopy(adata)
# ORI Works when loading from disk
scvi.model.SCVI.prepare_query_anndata(adata, "dataloader_model2")
# ORI This one still needs to be fixed.
scvi.model.SCVI.prepare_query_anndata(adata, model_census2)

# ORI Should work when setting up the AnnData correctly. scANVI with DataModule is not yet
# supported as DataModule can't take a labels_key.
scanvae = scvi.model.SCANVI.from_scvi_model(
    model_census2,
    adata=adata,
    unlabeled_category="Unknown",
    labels_key="cell_type",
)

# ORI - check it should work with a model initialized with AnnData. See below not fully working yet
model_census3 = scvi.model.SCVI.load("dataloader_model2", adata=adata)

scvi.model.SCVI.prepare_query_anndata(adata, "dataloader_model2")
query_model = scvi.model.SCVI.load_query_data(adata, "dataloader_model2")

scvi.model.SCVI.prepare_query_anndata(adata, model_census3)
query_model = scvi.model.SCVI.load_query_data(adata, model_census3)

# with open("model.pt", "rb") as f:
#     torch_model = torch.load(f)
#
#     adict = torch_model["attr_dict"]
#     params = adict["init_params_"]["non_kwargs"]
#
#     n_batch = adict["n_batch"]
#     n_extra_categorical_covs = adict["n_extra_categorical_covs"]
#     n_extra_continuous_covs = adict["n_extra_continuous_covs"]
#     n_labels = adict["n_labels"]
#     n_vars = adict["n_vars"]
#
#     latent_distribution = params["latent_distribution"]
#     dispersion = params["dispersion"]
#     n_hidden = params["n_hidden"]
#     dropout_rate = params["dropout_rate"]
#     gene_likelihood = params["gene_likelihood"]
#
#     model = scvi.model.SCVI(
#         n_layers=params["n_layers"],
#         n_latent=params["n_latent"],
#         gene_likelihood=params["gene_likelihood"],
#         encode_covariates=False,
#     )
#
#     module = model._module_cls(
#         n_input=n_vars,
#         n_batch=n_batch,
#         n_labels=n_labels,
#         n_continuous_cov=n_extra_continuous_covs,
#         n_cats_per_cov=None,
#         n_hidden=n_hidden,
#         n_latent=n_latent,
#         n_layers=n_layers,
#         dropout_rate=dropout_rate,
#         dispersion=dispersion,
#         gene_likelihood=gene_likelihood,
#         latent_distribution=latent_distribution,
#     )
#     model.module = module
#
#     model.module.load_state_dict(torch_model["model_state_dict"])
#
#     device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#
#     model.to_device(device)
#     model.module.eval()
#     model.is_trained = True

# We will now generate the cell embeddings for this model, using the `get_latent_representation`
# function available in scvi-tools.
# We can use another instance of the `ExperimentDataPipe` for the forward pass, so we don't need
# to load the whole dataset in memory.

# # Needs to have shuffle=False for inference
# datamodule_inference = CensusSCVIDataModule(
#     census["census_data"][experiment_name],
#     measurement_name="RNA",
#     X_name="raw",
#     obs_query=soma.AxisQuery(value_filter=obs_value_filter),
#     var_query=soma.AxisQuery(coords=(list(hv_idx),)),
#     batch_size=1024,
#     shuffle=False,
#     soma_chunk_size=50_000,
#     batch_keys=["dataset_id", "assay", "suspension_type", "donor_id"],
#     dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
# )
#
# # We can simply feed the datapipe to `get_latent_representation` to obtain the embeddings -
# # will take a while
# datapipe = datamodule_inference.datapipe
# dataloader = experiment_dataloader(datapipe, num_workers=0, persistent_workers=False)
# mapped_dataloader = (
#     datamodule_inference.on_before_batch_transfer(tensor, None) for tensor in dataloader
# )
# latent = model.get_latent_representation(dataloader=mapped_dataloader)
# emb_idx = datapipe._obs_joinids
#
# # We will now take a look at the UMAP for the generated embedding
# # (will be later comapred to what we got)
# adata = cellxgene_census.get_anndata(
#     census,
#     organism=experiment_name,
#     obs_value_filter=obs_value_filter,
# )
# obs_soma_joinids = adata.obs["soma_joinid"]
# obs_indexer = pd.Index(emb_idx)
# idx = obs_indexer.get_indexer(obs_soma_joinids)
# # Reindexing is necessary to ensure that the cells in the embedding match the
# # ones in the anndata object.
# adata.obsm["scvi"] = latent[idx]
#
# # Plot UMAP and save the figure for later check
# sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi")
# sc.tl.umap(adata, neighbors_key="scvi")
# sc.pl.umap(adata, color="dataset_id", title="SCVI")
#
#
# # Now return and add all the registry stuff that we will need
#
#
# # Now add the missing stuff from the current CZI implemenation in order for us to have the exact
# # same steps like the original way (except than setup_anndata)
