from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
import torch
from lightning.pytorch import LightningDataModule
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader

import scvi
from scvi.utils import dependencies

if TYPE_CHECKING:
    from typing import Any

    import numpy as np
    import pandas as pd


@dependencies("lamindb")
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
                        "unlabeled_category": self._unlabeled_category,
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


@dependencies("tiledbsoma")
@dependencies("tiledbsoma_ml")
class CensusSCVIDataModule(LightningDataModule):
    import tiledbsoma as soma

    """PyTorch Lightning DataModule for training scVI models from SOMA data

    Wraps a `tiledbsoma_ml.ExperimentDataset` to stream the results of a SOMA
    `ExperimentAxisQuery`, exposing a `DataLoader` to generate tensors ready for scVI model
    training. Also handles deriving the scVI batch label as a tuple of obs columns.
    """

    def __init__(
        self,
        query: soma.ExperimentAxisQuery,
        *args,
        batch_column_names: list[str] | None = None,
        batch_labels: list[str] | None = None,
        dataloader_kwargs: dict[str, Any] | None = None,
        **kwargs,
    ):
        """
        Args:

        query: tiledbsoma.ExperimentAxisQuery
                        Defines the desired result set from a SOMA Expeirement.
        *args, **kwargs:
        Additional arguments passed through to `tiledbsoma_ml.ExperimentDataset`.

        batch_column_names: List[str], optional
        List of obs column names, the tuple of which defines the scVI batch label
        (not to to be confused with a batch of training data). Defaults to
        `["dataset_id", "assay", "suspension_type", "donor_id"]`.

        batch_labels: List[str], optional
        List of possible values of the batch label, for mapping to label tensors. By default,
        this will be derived from the unique labels in the given query results (given
        `batch_column_names`), making the label mapping depend on the query. The `batch_labels`
        attribute in the `SCVIDataModule` used for training may be saved and here restored in
        another instance for a different query. That ensures the label mapping will be correct
        for the trained model, even if the second query doesn't return examples of every
        training batch label.

        dataloader_kwargs: dict, optional
        Keyword arguments passed to `tiledbsoma_ml.experiment_dataloader()`, e.g. `num_workers`.
        """
        super().__init__()
        self.query = query
        self.dataset_args = args
        self.dataset_kwargs = kwargs
        self.dataloader_kwargs = dataloader_kwargs if dataloader_kwargs is not None else {}
        self.batch_column_names = (
            batch_column_names
            if batch_column_names is not None
            else ["dataset_id", "assay", "suspension_type", "donor_id"]
        )
        self.batch_colsep = "//"
        self.batch_colname = "scvi_batch"
        # prepare LabelEncoder for the scVI batch label:
        #   1. read obs DataFrame for the whole query result set
        #   2. add scvi_batch column
        #   3. fit LabelEncoder to the scvi_batch column's unique values
        if batch_labels is None:
            obs_df = self.query.obs(column_names=self.batch_column_names).concat().to_pandas()
            self._add_batch_col(obs_df, inplace=True)
            batch_labels = obs_df[self.batch_colname].unique()
        self.batch_labels = batch_labels
        self.batch_encoder = LabelEncoder().fit(self.batch_labels)

    def setup(self, stage: str | None = None) -> None:
        # Instantiate the ExperimentDataset with the provided args and kwargs.
        import tiledbsoma_ml

        self.train_dataset = tiledbsoma_ml.ExperimentDataset(
            self.query,
            *self.dataset_args,
            obs_column_names=self.batch_column_names,
            **self.dataset_kwargs,
        )

    def train_dataloader(self) -> DataLoader:
        import tiledbsoma_ml

        return tiledbsoma_ml.experiment_dataloader(
            self.train_dataset,
            **self.dataloader_kwargs,
        )

    def _add_batch_col(self, obs_df: pd.DataFrame, inplace: bool = False):
        # synthesize a new column for obs_df by concatenating the self.batch_column_names columns
        if not inplace:
            obs_df = obs_df.copy()
        obs_df[self.batch_colname] = (
            obs_df[self.batch_column_names].astype(str).agg(self.batch_colsep.join, axis=1)
        )
        return obs_df

    def on_before_batch_transfer(
        self,
        batch,
        dataloader_idx: int,
    ) -> dict[str, torch.Tensor | None]:
        # DataModule hook: transform the ExperimentDataset data batch
        # (X: ndarray, obs_df: DataFrame)
        # into X & batch variable tensors for scVI (using batch_encoder on scvi_batch)
        batch_X, batch_obs = batch
        self._add_batch_col(batch_obs, inplace=True)
        return {
            "X": torch.from_numpy(batch_X).float(),
            "batch": torch.from_numpy(
                self.batch_encoder.transform(batch_obs[self.batch_colname])
            ).unsqueeze(1),
            "labels": torch.empty(0),
        }

    # scVI code expects these properties on the DataModule:

    @property
    def n_obs(self) -> int:
        return len(self.query.obs_joinids())

    @property
    def n_vars(self) -> int:
        return len(self.query.var_joinids())

    @property
    def n_batch(self) -> int:
        return len(self.batch_encoder.classes_)
