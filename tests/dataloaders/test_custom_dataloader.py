from __future__ import annotations

from pprint import pprint
from time import time
from typing import TYPE_CHECKING

import psutil
from lightning.pytorch import LightningDataModule
from torch.utils.data import DataLoader

import scvi
from scvi.utils import dependencies

if TYPE_CHECKING:
    import lamindb as ln
    import numpy as np


# @pytest.mark.custom_dataloader
@dependencies("lamindb")
def test_lamindb_dataloader_scvi_scanvi(save_path: str):
    import lamindb as ln

    class MappedCollectionDataModule(LightningDataModule):
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
                    "size_factor": {
                        "data_registry": {},
                        "state_registry": {},
                        "summary_stats": {},
                    },
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
                BATCH_KEY: batch[self._batch_key][:, None]
                if self._batch_key is not None
                else None,
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
