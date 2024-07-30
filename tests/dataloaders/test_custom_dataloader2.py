from __future__ import annotations

import os

import cellxgene_census
import pandas as pd
import scanpy as sc
import tiledbsoma as soma
import torch
from cellxgene_census.experimental.pp import highly_variable_genes

import scvi
from scvi.dataloaders._custom_dataloader import CensusSCVIDataModule, experiment_dataloader
from scvi.model import SCVI

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
# This is a new func to implement
SCVI.setup_datamodule(datamodule)
#
model = SCVI(n_layers=n_layers, n_latent=n_latent, gene_likelihood="nb", encode_covariates=False)


# The CZI data module is a refined data module while SCVI is a lighting datamodule
# Altough this is only 1 epoch it will take few mins on local machine
model.train(
    datamodule=datamodule,
    max_epochs=max_epochs,
    batch_size=batch_size,
    train_size=train_size,
    early_stopping=False,
)

# We can now save the trained model. As of the current writing date (June 2024),
# scvi-tools doesn't support saving a model that wasn't generated through an AnnData loader,
# so we'll use some custom code:
model_state_dict = model.module.state_dict()
var_names = hv_idx.to_numpy()
user_attributes = model._get_user_attributes()
user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

user_attributes.update(
    {
        "n_batch": datamodule.n_batch,
        "n_extra_categorical_covs": 0,
        "n_extra_continuous_covs": 0,
        "n_labels": 1,
        "n_vars": datamodule.n_vars,
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

# Saving the model the original way
save_dir = "/Users/orikr/runs/290724/"  # tempfile.TemporaryDirectory()
model_dir = os.path.join(save_dir, "scvi_czi_model")
model.save(model_dir, overwrite=True)


# We will now load the model back and use it to generate cell embeddings (the latent space),
# which can then be used for further analysis. Note that we still need to use some custom code for
# loading the model, which includes loading the parameters from the `attr_dict` node stored in
# the model.
with open("model.pt", "rb") as f:
    torch_model = torch.load(f)

    adict = torch_model["attr_dict"]
    params = adict["init_params_"]["non_kwargs"]

    n_batch = adict["n_batch"]
    n_extra_categorical_covs = adict["n_extra_categorical_covs"]
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
        n_layers=n_layers,
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

# We will now generate the cell embeddings for this model, using the `get_latent_representation`
# function available in scvi-tools.
# We can use another instance of the `ExperimentDataPipe` for the forward pass, so we don't need
# to load the whole dataset in memory.

# Needs to have shuffle=False for inference
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

# We can simply feed the datapipe to `get_latent_representation` to obtain the embeddings -
# will take a while
datapipe = datamodule_inference.datapipe
dataloader = experiment_dataloader(datapipe, num_workers=0, persistent_workers=False)
mapped_dataloader = (
    datamodule_inference.on_before_batch_transfer(tensor, None) for tensor in dataloader
)
latent = model.get_latent_representation(dataloader=mapped_dataloader)
emb_idx = datapipe._obs_joinids

# We will now take a look at the UMAP for the generated embedding
# (will be later comapred to what we got)
adata = cellxgene_census.get_anndata(
    census,
    organism=experiment_name,
    obs_value_filter=obs_value_filter,
)
obs_soma_joinids = adata.obs["soma_joinid"]
obs_indexer = pd.Index(emb_idx)
idx = obs_indexer.get_indexer(obs_soma_joinids)
# Reindexing is necessary to ensure that the cells in the embedding match the
# ones in the anndata object.
adata.obsm["scvi"] = latent[idx]

# Plot UMAP and save the figure for later check
sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi")
sc.tl.umap(adata, neighbors_key="scvi")
sc.pl.umap(adata, color="dataset_id", title="SCVI")


# Now return and add all the registry stuff that we will need


# Now add the missing stuff from the current CZI implemenation in order for us to have the exact
# same steps like the original way (except than setup_anndata)
