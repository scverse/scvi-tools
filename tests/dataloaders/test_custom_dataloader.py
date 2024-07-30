from __future__ import annotations

import os

import numpy as np
import scanpy as sc

import scvi
from scvi.data import _constants, synthetic_iid
from scvi.model import SCVI

# We will now create the SCVI model object:
# Its parameters:
n_layers = 1
n_latent = 10
batch_size = 1024
train_size = 0.9
max_epochs = 1


# COMAPRE TO THE ORIGINAL METHOD!!! - use the same data!!!
# We first create a registry using the orignal way of anndata in order to compare and add
# what is missing
adata = synthetic_iid()
adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
SCVI.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="labels",
    size_factor_key="size_factor",
)
#
model_orig = SCVI(adata, n_latent=n_latent)
model_orig.train(1, check_val_every_n_epoch=1, train_size=0.5)

# Saving the model
save_dir = "/Users/orikr/runs/290724/"  # tempfile.TemporaryDirectory()
model_dir = os.path.join(save_dir, "scvi_orig_model")
model_orig.save(model_dir, overwrite=True)

# Loading the model (just as a compariosn)
model_orig_loaded = scvi.model.SCVI.load(model_dir, adata=adata)

# Obtaining model outputs
SCVI_LATENT_KEY = "X_scVI"
latent = model_orig.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
# latent.shape

# You can see all necessary entries and the structure at
adata_manager = model_orig.adata_manager
model_orig.view_anndata_setup(hide_state_registries=True)
# adata_manager.get_state_registry(SCVI.REGISTRY_KEYS.X_KEY).to_dict()
adata_manager.registry[_constants._FIELD_REGISTRIES_KEY]

# Plot UMAP and save the figure for later check
sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi")
sc.tl.umap(adata, neighbors_key="scvi")
sc.pl.umap(adata, color="dataset_id", title="SCVI")

# Now return and add all the registry stuff that we will need

# Now add the missing stuff from the current CZI implemenation in order for us to have the exact
# same steps like the original way (except than setup_anndata)
