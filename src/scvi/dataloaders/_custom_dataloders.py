from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
from lightning.pytorch import LightningDataModule
from torch.utils.data import DataLoader

import scvi

if TYPE_CHECKING:
    import lamindb as ln
    import numpy as np


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
