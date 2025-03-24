from __future__ import annotations

import os
from pprint import pprint
from time import time

import pytest

import scvi
from scvi.dataloaders import MappedCollectionDataModule
from scvi.utils import dependencies


@pytest.mark.dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_scanvi(save_path: str):
    os.system("lamin init --storage ./test-registries")
    import lamindb as ln

    ln.setup.init(name="lamindb_instance_name", storage=save_path)

    # a test for mapped collection
    collection = ln.Collection.get(name="covid_normal_lung")
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
