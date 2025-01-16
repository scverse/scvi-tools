from __future__ import annotations

import os
import sys
from pprint import pprint
from time import time

import numpy as np
import pandas as pd
import psutil
import pytest
import scanpy as sc
import tqdm
from lightning.pytorch import LightningDataModule
from torch.utils.data import DataLoader

import scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.utils import dependencies


class MappedCollectionDataModule(LightningDataModule):
    import lamindb as ln

    def __init__(
        self,
        collection: ln.Collection,
        batch_key: str | None = None,
        label_key: str | None = None,
        batch_size: int = 128,
        **kwargs,
    ):
        self._batch_size = batch_size
        self._batch_key = batch_key
        self._label_key = label_key
        self._parallel = kwargs.pop("parallel", True)
        # here we initialize MappedCollection to use in a pytorch DataLoader
        self._dataset = collection.mapped(
            obs_keys=self._batch_key, parallel=self._parallel, **kwargs
        )
        # need by scvi and lightning.pytorch
        self._log_hyperparams = False
        self.allow_zero_length_dataloader_with_multiple_devices = False

    def close(self):
        self._dataset.close()

    def setup(self, stage):
        pass

    def train_dataloader(self):
        return self._create_dataloader(shuffle=True)

    def inference_dataloader(self):
        """Dataloader for inference with `on_before_batch_transfer` applied."""
        dataloader = self._create_dataloader(shuffle=False, batch_size=4096)
        return self._InferenceDataloader(dataloader, self.on_before_batch_transfer)

    def _create_dataloader(self, shuffle, batch_size=None):
        if self._parallel:
            num_workers = psutil.cpu_count() - 1
            worker_init_fn = self._dataset.torch_worker_init_fn
        else:
            num_workers = 0
            worker_init_fn = None
        if batch_size is None:
            batch_size = self._batch_size
        return DataLoader(
            self._dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            num_workers=num_workers,
            worker_init_fn=worker_init_fn,
        )

    @property
    def n_obs(self) -> int:
        return self._dataset.n_obs

    @property
    def var_names(self) -> int:
        return self._dataset.var_joint

    @property
    def n_vars(self) -> int:
        return self._dataset.n_vars

    @property
    def n_batch(self) -> int:
        if self._batch_key is None:
            return 1
        return len(self._dataset.encoders[self._batch_key])

    @property
    def n_labels(self) -> int:
        if self._label_key is None:
            return 1
        return len(self._dataset.encoders[self._label_key])

    @property
    def labels(self) -> np.ndarray:
        return self._dataset[self._label_key]

    @property
    def registry(self) -> dict:
        return {
            "scvi_version": scvi.__version__,
            "model_name": "SCVI",
            "setup_args": {
                "layer": None,
                "batch_key": self._batch_key,
                "labels_key": self._label_key,
                "size_factor_key": None,
                "categorical_covariate_keys": None,
                "continuous_covariate_keys": None,
            },
            "field_registries": {
                "X": {
                    "data_registry": {"attr_name": "X", "attr_key": None},
                    "state_registry": {
                        "n_obs": self.n_obs,
                        "n_vars": self.n_vars,
                        "column_names": self.var_names,
                    },
                    "summary_stats": {"n_vars": self.n_vars, "n_cells": self.n_obs},
                },
                "batch": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_batch"},
                    "state_registry": {
                        "categorical_mapping": self.batch_keys,
                        "original_key": self._batch_key,
                    },
                    "summary_stats": {"n_batch": self.n_batch},
                },
                "labels": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_labels"},
                    "state_registry": {
                        "categorical_mapping": self.label_keys,
                        "original_key": self._label_key,
                        "unlabeled_category": "unlabeled",
                    },
                    "summary_stats": {"n_labels": self.n_labels},
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
            "setup_method_name": "setup_anndata",
        }

    @property
    def batch_keys(self) -> int:
        if self._batch_key is None:
            return None
        return self._dataset.encoders[self._batch_key]

    @property
    def label_keys(self) -> int:
        if self._label_key is None:
            return None
        return self._dataset.encoders[self._label_key]

    def on_before_batch_transfer(
        self,
        batch,
        dataloader_idx,
    ):
        X_KEY: str = "X"
        BATCH_KEY: str = "batch"
        LABEL_KEY: str = "labels"

        return {
            X_KEY: batch["X"].float(),
            BATCH_KEY: batch[self._batch_key][:, None] if self._batch_key is not None else None,
            LABEL_KEY: 0,
        }

    class _InferenceDataloader:
        """Wrapper to apply `on_before_batch_transfer` during iteration."""

        def __init__(self, dataloader, transform_fn):
            self.dataloader = dataloader
            self.transform_fn = transform_fn

        def __iter__(self):
            for batch in self.dataloader:
                yield self.transform_fn(batch, dataloader_idx=None)

        def __len__(self):
            return len(self.dataloader)


@pytest.mark.custom_dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_scanvi(save_path: str = "."):
    import lamindb as ln

    # a test for mapped collection
    collection = ln.Collection.get(name="covid_normal_lung")
    datamodule = MappedCollectionDataModule(
        collection, batch_key="assay", batch_size=128, join="inner"
    )
    model = scvi.model.SCVI(adata=None, registry=datamodule.registry)
    pprint(model.summary_stats)
    pprint(model.module)
    dataloader = datamodule.inference_dataloader()
    adata = collection.load()  # try to compare this in regular settings

    # Using regular adata laoder
    # setup large
    SCVI.setup_anndata(
        adata,
        batch_key="assay",
    )
    model_reg = SCVI(adata)
    start_time = time()
    model_reg.train(max_epochs=10, batch_size=1024)
    time_reg = time() - start_time
    print(time_reg)

    start_time = time()
    model.train(max_epochs=10, batch_size=1024, datamodule=datamodule)
    time_lamin = time() - start_time
    print(time_lamin)

    _ = model.get_elbo(dataloader=dataloader)
    _ = model.get_marginal_ll(dataloader=dataloader)
    _ = model.get_reconstruction_error(dataloader=dataloader)
    _ = model.get_latent_representation(dataloader=dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True)
    model_query = model.load_query_data(
        adata=False, reference_model="lamin_model", registry=datamodule.registry
    )
    model_query.train(max_epochs=1, datamodule=datamodule)
    _ = model_query.get_elbo(dataloader=dataloader)
    _ = model_query.get_marginal_ll(dataloader=dataloader)
    _ = model_query.get_reconstruction_error(dataloader=dataloader)
    _ = model_query.get_latent_representation(dataloader=dataloader)

    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata=adata, reference_model="lamin_model")
    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata)
    model_query_adata.train(max_epochs=1)
    _ = model_query_adata.get_elbo()
    _ = model_query_adata.get_marginal_ll()
    _ = model_query_adata.get_reconstruction_error()
    _ = model_query_adata.get_latent_representation()
    _ = model_query_adata.get_latent_representation(dataloader=dataloader)

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


@pytest.mark.custom_dataloader
def test_czi_custom_dataloader_scvi(save_path: str = "."):
    # should be ready for importing the cloned branch on a remote machine that runs github action
    sys.path.insert(
        0,
        "/home/runner/work/scvi-tools/scvi-tools/cellxgene-census/api/python/cellxgene_census/src",
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

    model = scvi.model.SCVI(adata_orig)
    model.train(max_epochs=1)

    dataloader = model._make_data_loader(adata_orig)
    _ = model.get_elbo(dataloader=dataloader)
    _ = model.get_marginal_ll(dataloader=dataloader)
    _ = model.get_reconstruction_error(dataloader=dataloader)
    _ = model.get_latent_representation(dataloader=dataloader)

    scvi.model.SCVI.prepare_query_anndata(adata_orig, reference_model=model)
    scvi.model.SCVI.load_query_data(adata_orig, reference_model=model)

    user_attributes = model._get_user_attributes()
    pprint(user_attributes)

    scvi.model._scvi.SCVI.setup_datamodule(datamodule)  # takes time
    model_census = scvi.model.SCVI(
        registry=datamodule.registry,
        gene_likelihood="nb",
        encode_covariates=False,
    )

    pprint(datamodule.registry)

    max_epochs = 1

    model_census.train(
        datamodule=datamodule,
        max_epochs=max_epochs,
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
def test_czi_custom_dataloader_scanvi(save_path: str = "."):
    # should be ready for importing the cloned branch on a remote machine that runs github action
    sys.path.insert(
        0,
        "/home/runner/work/scvi-tools/scvi-tools/cellxgene-census/api/python/cellxgene_census/src",
    )
    sys.path.insert(0, "src")
    import cellxgene_census
    import tiledbsoma as soma
    from cellxgene_census.experimental.ml.datamodule import CensusSCVIDataModule
    # from cellxgene_census.experimental.pp import highly_variable_genes

    # this test checks the local custom dataloder made by CZI and run several tests with it
    census = cellxgene_census.open_soma(census_version="stable")

    experiment_name = "mus_musculus"
    obs_value_filter = (
        'is_primary_data == True and tissue_general in ["kidney","liver"] and nnz >= 3000'
    )

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
    label_keys = ["tissue_general"]
    datamodule = CensusSCVIDataModule(
        census["census_data"][experiment_name],
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(list(hv_idx),)),
        batch_size=1024,
        shuffle=True,
        batch_keys=batch_keys,
        label_keys=label_keys,
        unlabeled_category="label_0",
        dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
    )

    # table of genes should be filtered by soma_joinid - but we should keep the encoded indexes
    # This is nice to have and might be uses in the downstream analysis
    # var_df = census["census_data"][experiment_name].ms["RNA"].var.read().concat().to_pandas()
    # var_df = var_df.loc[var_df.soma_joinid.isin(
    #     list(datamodule.datapipe.var_query.coords[0] if datamodule.datapipe.var_query is not None
    #          else range(datamodule.n_vars)))]

    # scvi.model._scvi.SCVI.setup_datamodule(datamodule)  # takes time
    # model_census = scvi.model.SCVI(
    #     registry=datamodule.registry,
    #     gene_likelihood="nb",
    #     encode_covariates=False,
    # )
    #
    # pprint(datamodule.registry)
    #
    max_epochs = 1
    #
    # model_census.train(
    #     datamodule=datamodule,
    #     max_epochs=max_epochs,
    #     early_stopping=False,
    # )

    scvi.model.SCANVI.setup_datamodule(datamodule)
    pprint(datamodule.registry)
    model = scvi.model.SCANVI(registry=datamodule.registry, datamodule=datamodule)
    model.view_anndata_setup(datamodule)
    adata_manager = model.adata_manager
    pprint(adata_manager.registry)
    model.train(
        datamodule=datamodule, max_epochs=max_epochs, train_size=0.5, check_val_every_n_epoch=1
    )
    # logged_keys = model.history.keys()
    # assert len(model._labeled_indices) == sum(adata.obs["labels"] != "label_0")
    # assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == "label_0")
    # assert "elbo_validation" in logged_keys
    # assert "reconstruction_loss_validation" in logged_keys
    # assert "kl_local_validation" in logged_keys
    # assert "elbo_train" in logged_keys
    # assert "reconstruction_loss_train" in logged_keys
    # assert "kl_local_train" in logged_keys
    # assert "validation_classification_loss" in logged_keys
    # assert "validation_accuracy" in logged_keys
    # assert "validation_f1_score" in logged_keys
    # assert "validation_calibration_error" in logged_keys
    # adata2 = synthetic_iid()
    # predictions = model.predict(adata2, indices=[1, 2, 3])
    # assert len(predictions) == 3
    # model.predict()
    # df = model.predict(adata2, soft=True)
    # assert isinstance(df, pd.DataFrame)
    # model.predict(adata2, soft=True, indices=[1, 2, 3])
    # model.get_normalized_expression(adata2)
    # model.differential_expression(groupby="labels", group1="label_1")
    # model.differential_expression(groupby="labels", group1="label_1", group2="label_2")


@pytest.mark.custom_dataloader
@dependencies("lamindb")
def test_scdataloader_custom_dataloader_scvi(save_path: str = "."):
    import lamindb as ln

    os.system("lamin init --storage ~/scdataloader2 --schema bionty")
    from scdataloader import Collator, DataModule, SimpleAnnDataset

    # import bionty as bt
    # from scdataloader import utils
    from scdataloader.preprocess import (
        LaminPreprocessor,
        additional_postprocess,
        additional_preprocess,
    )

    # import tiledbsoma as soma
    from scdataloader.utils import populate_my_ontology
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
    max_epochs = 1

    # Try the lamindb dataloder on a trained scvi-model with adata
    # adata = adata.copy()
    scvi.model.SCVI.setup_anndata(adata_orig, batch_key="cell_type_ontology_term_id")
    model = scvi.model.SCVI(adata_orig)
    model.train(max_epochs=max_epochs)
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
    #     gene_likelihood="nb",
    #     encode_covariates=False,
    # )
    #
    # pprint(datamodule.registry)
    #
    # model_lamindb.train(
    #     datamodule=datamodule,
    #     max_epochs=max_epochs,
    #     early_stopping=False,
    # )
    # We have to create a registry without setup_anndata that contains the same elements
    # The other way will be to fill the model ,LIKE IN CELLXGENE NOTEBOOK
    # need to pass here new object of registry taht contains everything we will need
