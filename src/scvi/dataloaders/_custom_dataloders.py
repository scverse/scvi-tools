from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import torch
from lightning.pytorch import LightningDataModule
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader

import scvi
from scvi.model._utils import parse_device_args
from scvi.utils import dependencies

if TYPE_CHECKING:
    from typing import Any

    import lamindb as ln
    import pandas as pd
    import tiledbsoma as soma


class MappedCollectionDataModule(LightningDataModule):
    @dependencies("lamindb")
    def __init__(
        self,
        collection: ln.Collection,
        batch_key: str | None = None,
        label_key: str | None = None,
        unlabeled_category: str | None = "Unknown",
        batch_size: int = 128,
        collection_val: ln.Collection | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        shuffle: bool = True,
        model_name: str = "SCVI",
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        super().__init__()
        self._batch_size = batch_size
        self._batch_key = batch_key
        self._label_key = label_key
        self.model_name = model_name
        self.shuffle = shuffle
        self.unlabeled_category = unlabeled_category
        self._parallel = kwargs.pop("parallel", True)
        self.labels_ = None
        self._categorical_covariate_keys = categorical_covariate_keys
        self._continuous_covariate_keys = continuous_covariate_keys

        # here we initialize MappedCollection to use in a pytorch DataLoader
        obs_keys = self._batch_key  # we must have batch keys
        if self._label_key is not None:
            obs_keys = [obs_keys] + [self._label_key]
        if self._categorical_covariate_keys is not None:
            obs_keys = (
                obs_keys + self._categorical_covariate_keys
                if type(obs_keys).__name__ == "list"
                else [obs_keys] + self._categorical_covariate_keys
            )
        if self._continuous_covariate_keys is not None:
            obs_keys = (
                obs_keys + self._continuous_covariate_keys
                if type(obs_keys).__name__ == "list"
                else [obs_keys] + self._continuous_covariate_keys
            )

        self._dataset = collection.mapped(obs_keys=obs_keys, parallel=self._parallel, **kwargs)
        if collection_val is not None:
            self._validset = collection_val.mapped(
                obs_keys=obs_keys, parallel=self._parallel, **kwargs
            )
        else:
            self._validset = None

        # generate encodings
        if self._label_key is not None:
            self.labels_ = self._dataset.get_merged_labels(self._label_key).astype(str)
        if self._categorical_covariate_keys is not None:
            self.categorical_covariate_keys_ = [
                self._dataset.encoders[cat_cov_key]
                for cat_cov_key in self._categorical_covariate_keys
            ]

        # need by scvi and lightning.pytorch
        self._log_hyperparams = False
        self.allow_zero_length_dataloader_with_multiple_devices = False
        _, _, self.device = parse_device_args(
            accelerator=accelerator, devices=device, return_device="torch"
        )

    def close(self):
        self._dataset.close()
        self._validset.close()

    def train_dataloader(self) -> DataLoader:
        return self._create_dataloader(shuffle=self.shuffle)

    def val_dataloader(self) -> DataLoader:
        return self._create_dataloader_val(shuffle=self.shuffle)

    def inference_dataloader(self, shuffle=False, batch_size=4096, indices=None):
        """Dataloader for inference with `on_before_batch_transfer` applied."""
        dataloader = self._create_dataloader(shuffle, batch_size, indices)
        return self._InferenceDataloader(dataloader, self.on_before_batch_transfer)

    def _create_dataloader(self, shuffle, batch_size=None, indices=None):
        if self._parallel:
            num_workers = os.cpu_count() - 1
            worker_init_fn = self._dataset.torch_worker_init_fn
        else:
            num_workers = 0
            worker_init_fn = None
        if batch_size is None:
            batch_size = self._batch_size
        if indices is not None:
            dataset = self._dataset[indices]  # TODO find a better way
        else:
            dataset = self._dataset
        return DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            num_workers=num_workers,
            worker_init_fn=worker_init_fn,
        )

    def _create_dataloader_val(self, shuffle, batch_size=None, indices=None):
        if self._validset is not None:
            if self._parallel:
                num_workers = os.cpu_count() - 1
                worker_init_fn = self._validset.torch_worker_init_fn
            else:
                num_workers = 0
                worker_init_fn = None
            if batch_size is None:
                batch_size = self._batch_size
            if indices is not None:
                validset = self._validset[indices]  # TODO find a better way
            else:
                validset = self._validset
            return DataLoader(
                validset,
                batch_size=batch_size,
                shuffle=shuffle,
                num_workers=num_workers,
                worker_init_fn=worker_init_fn,
            )
        else:
            pass

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
            return 0
        return len(self.labels)

    @property
    def labels(self) -> np.ndarray:
        if self._label_key is None:
            return None
        combined = np.concatenate(
            (list(self._dataset.encoders[self._label_key].keys()), [self.unlabeled_category])
        )
        unique_values, idx = np.unique(combined, return_index=True)
        unique_values = unique_values[np.argsort(idx)]
        return unique_values.astype(object)

    @property
    def unlabeled_category(self) -> str:
        """String assigned to unlabeled cells."""
        if not hasattr(self, "_unlabeled_category"):
            raise AttributeError("`unlabeled_category` not set.")
        return self._unlabeled_category

    @unlabeled_category.setter
    def unlabeled_category(self, value: str | None):
        if not (value is None or isinstance(value, str)):
            raise ValueError("`unlabeled_category` must be a string or None.")
        self._unlabeled_category = value

    @property
    def extra_categorical_covs(self) -> dict:
        if self._categorical_covariate_keys is None:
            out = {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_categorical_covs": 0},
            }
        else:
            # TODO need to adjust this mapping
            mapping = dict(
                zip(
                    self._categorical_covariate_keys,
                    [np.array(list(x.values())) for x in self.categorical_covariate_keys_],
                    strict=False,
                )
            )
            out = {
                "data_registry": {"attr_key": "_scvi_extra_categorical_covs", "attr_name": "obsm"},
                "state_registry": {
                    "field_keys": self._categorical_covariate_keys,
                    "mapping": mapping,
                    "n_cats_per_key": [len(mapping[map]) for map in mapping.keys()],
                },
                "summary_stats": {
                    "n_extra_categorical_covs": len(self._categorical_covariate_keys)
                },
            }
        return out

    @property
    def extra_continuous_covs(self) -> dict:
        if self._continuous_covariate_keys is None:
            out = {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_continuous_covs": 0},
            }
        else:
            out = {
                "data_registry": {"attr_key": "_scvi_extra_continuous_covs", "attr_name": "obsm"},
                "state_registry": {
                    "columns": np.array(self._continuous_covariate_keys, dtype=object)
                },
                "summary_stats": {"n_extra_continuous_covs": len(self._continuous_covariate_keys)},
            }
        return out

    @property
    def registry(self) -> dict:
        return {
            "scvi_version": scvi.__version__,
            "model_name": self.model_name,
            "setup_args": {
                "layer": None,
                "batch_key": self._batch_key,
                "labels_key": self._label_key,
                "size_factor_key": None,
                "categorical_covariate_keys": self._categorical_covariate_keys,
                "continuous_covariate_keys": self._continuous_covariate_keys,
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
                        "categorical_mapping": self.batch_labels,
                        "original_key": self._batch_key,
                    },
                    "summary_stats": {"n_batch": self.n_batch},
                },
                "labels": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_labels"},
                    "state_registry": {
                        "categorical_mapping": self.labels,
                        "original_key": self._label_key,
                        "unlabeled_category": self.unlabeled_category,
                    },
                    "summary_stats": {"n_labels": self.n_labels},
                },
                "size_factor": {
                    "data_registry": {},
                    "state_registry": {},
                    "summary_stats": {},
                },
                "extra_categorical_covs": self.extra_categorical_covs,
                "extra_continuous_covs": self.extra_continuous_covs,
            },
            "setup_method_name": "setup_datamodule",
        }

    @property
    def batch_labels(self) -> int | None:
        if self._batch_key is None:
            return None
        return self._dataset.encoders[self._batch_key]

    @property
    def label_keys(self) -> int | None:
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
        CAT_COVS_KEY: str = "extra_categorical_covs"
        CONT_COVS_KEY: str = "extra_continuous_covs"

        return {
            X_KEY: batch["X"].float(),
            BATCH_KEY: batch[self._batch_key][:, None] if self._batch_key is not None else None,
            LABEL_KEY: batch[self._label_key][:, None] if self._label_key is not None else 0,
            CAT_COVS_KEY: torch.cat(
                [batch[k][:, None] for k in self._categorical_covariate_keys], dim=1
            )
            if self._categorical_covariate_keys is not None
            else None,
            CONT_COVS_KEY: torch.cat(
                [batch[k][:, None] for k in self._continuous_covariate_keys], dim=1
            )
            if self._continuous_covariate_keys is not None
            else None,
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


class TileDBDataModule(LightningDataModule):
    """PyTorch Lightning DataModule for training scVI models from SOMA data

    Wraps a `tiledbsoma_ml.ExperimentDataset` to stream the results of a SOMA
    `ExperimentAxisQuery`, exposing a `DataLoader` to generate tensors ready for scVI model
    training. Also handles deriving the scVI batch label as a tuple of obs columns.
    """

    @dependencies("tiledbsoma")
    def __init__(
        self,
        query: soma.ExperimentAxisQuery,
        *args,
        batch_column_names: list[str] | None = None,
        batch_labels: list[str] | None = None,
        label_keys: list[str] | None = None,
        unlabeled_category: str | None = "Unknown",
        train_size: float | None = 1.0,
        split_seed: int | None = None,
        dataloader_kwargs: dict[str, Any] | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        model_name: str = "SCVI",
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """
        Args:

        query: tiledbsoma.ExperimentAxisQuery
                        Defines the desired result set from a SOMA Experiment.
        *args, **kwargs:
        Additional arguments passed through to `tiledbsoma_ml.ExperimentDataset`.

        batch_column_names: List[str], optional
        List of obs column names, the tuple of which defines the scVI batch label
        (not to be confused with a batch of training data).

        batch_labels: List[str], optional
        List of possible values of the batch label, for mapping to label tensors. By default,
        this will be derived from the unique labels in the given query results (given
        `batch_column_names`), making the label mapping depend on the query. The `batch_labels`
        attribute in the `TileDBDataModule` used for training may be saved and here restored in
        another instance for a different query. That ensures the label mapping will be correct
        for the trained model, even if the second query doesn't return examples of every
        training batch label.

        label_keys
            List of obs column names concatenated to form the label column.
        unlabeled_category
            Value used for unlabeled cells in `labels_key` used to set up CZI datamodule with scvi.

        train_size
            Fraction of data to use for training.
        split_seed
            Seed for data split.

        dataloader_kwargs: dict, optional
        %(param_accelerator)s
        %(param_device)s

        model_name
            The SCVI-Tools Model we are running
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s

        Keyword arguments passed to `tiledbsoma_ml.experiment_dataloader()`, e.g. `num_workers`.
        """
        super().__init__()
        self.query = query
        self.dataset_args = args
        self.dataset_kwargs = kwargs
        self.dataloader_kwargs = dataloader_kwargs if dataloader_kwargs is not None else {}
        self.train_size = train_size
        self.split_seed = split_seed
        self.model_name = model_name

        # deal with labels if needed
        self.unlabeled_category = unlabeled_category
        self.label_keys = label_keys
        self.labels_colsep = "//"
        self.label_colname = "_scvi_labels"
        self.labels = None
        self.label_encoder = None
        self.labels_ = None
        self._categorical_covariate_keys = categorical_covariate_keys
        self._continuous_covariate_keys = continuous_covariate_keys
        self.categ_cov_colsep = "//"
        self._categorical_covariate_colname = "_scvi_cat_cov"

        # deal with batches
        self.batch_column_names = batch_column_names
        self.batch_colsep = "//"
        self.batch_colname = "_scvi_batch"
        # prepare LabelEncoder for the scVI batch label:
        #   1. read obs DataFrame for the whole query result set
        #   2. add scvi_batch column
        #   3. fit LabelEncoder to the scvi_batch column's unique values
        if batch_labels is None:
            cols_sel = (
                self.batch_column_names
                if self.label_keys is None
                else self.batch_column_names + self.label_keys
            )
            cols_sel = (
                cols_sel
                if self._categorical_covariate_keys is None
                else cols_sel + self._categorical_covariate_keys
            )
            cols_sel = (
                cols_sel
                if self._continuous_covariate_keys is None
                else cols_sel + self._continuous_covariate_keys
            )

            obs_df = self.query.obs(column_names=cols_sel).concat().to_pandas()
            obs_df = obs_df[cols_sel]
            self._add_batch_col(obs_df, inplace=True)
            batch_labels = obs_df[self.batch_colname].unique()
        self.batch_labels = batch_labels
        self.batch_encoder = LabelEncoder().fit(self.batch_labels)

        if label_keys is not None:
            obs_label_df = self.query.obs(column_names=self.label_keys).concat().to_pandas()
            obs_label_df = obs_label_df[self.label_keys]
            self._add_label_col(obs_label_df, inplace=True)
            labels = obs_label_df[self.label_colname].unique()
            self.labels = labels
            self.label_encoder = LabelEncoder().fit(self.labels)
            self.labels_ = obs_label_df[self.label_colname].values

        if categorical_covariate_keys is not None:
            obs_categ_cov_df = (
                self.query.obs(column_names=self._categorical_covariate_keys).concat().to_pandas()
            )
            obs_categ_cov_df = obs_categ_cov_df[self._categorical_covariate_keys]
            self._add_categ_cov_col(obs_categ_cov_df, inplace=True)
            categ_cov = obs_categ_cov_df[self._categorical_covariate_colname].unique()
            self.categ_cov = categ_cov
            self.categ_cov_encoder = LabelEncoder().fit(self.categ_cov)
            self.categ_cov_ = obs_categ_cov_df[self._categorical_covariate_colname].values

        _, _, self.device = parse_device_args(
            accelerator=accelerator, devices=device, return_device="torch"
        )

    @dependencies("tiledbsoma_ml")
    def setup(self, stage: str | None = None) -> None:
        # Instantiate the ExperimentDataset with the provided args and kwargs.
        from tiledbsoma_ml import ExperimentDataset

        cols_sel = (
            self.batch_column_names
            if self.label_keys is None
            else self.batch_column_names + self.label_keys
        )
        cols_sel = (
            cols_sel
            if self._categorical_covariate_keys is None
            else cols_sel + self._categorical_covariate_keys
        )
        cols_sel = (
            cols_sel
            if self._continuous_covariate_keys is None
            else cols_sel + self._continuous_covariate_keys
        )

        self.train_dataset = ExperimentDataset(
            self.query,
            *self.dataset_args,
            obs_column_names=cols_sel,
            **self.dataset_kwargs,
        )

        if self.validation_size > 0.0:
            datapipes = self.train_dataset.random_split(
                self.train_size, self.validation_size, seed=self.split_seed
            )
            self.train_dataset = datapipes[0]
            self.val_dataset = datapipes[1]
        else:
            self.val_dataset = None

    @dependencies("tiledbsoma_ml")
    def train_dataloader(self) -> DataLoader:
        from tiledbsoma_ml import experiment_dataloader

        return experiment_dataloader(
            self.train_dataset,
            **self.dataloader_kwargs,
        )

    @dependencies("tiledbsoma_ml")
    def val_dataloader(self) -> DataLoader:
        from tiledbsoma_ml import experiment_dataloader

        if self.val_dataset is not None:
            return experiment_dataloader(
                self.val_dataset,
                **self.dataloader_kwargs,
            )
        else:
            pass

    def _add_batch_col(self, obs_df: pd.DataFrame, inplace: bool = False):
        # synthesize a new column for obs_df by concatenating the self.batch_column_names columns
        if not inplace:
            obs_df = obs_df.copy()
        obs_df[self.batch_colname] = (
            obs_df[self.batch_column_names].astype(str).agg(self.batch_colsep.join, axis=1)
        )
        if self.labels is not None:
            obs_df[self.label_colname] = (
                obs_df[self.label_keys].astype(str).agg(self.labels_colsep.join, axis=1)
            )
        return obs_df

    def _add_label_col(self, obs_label_df: pd.DataFrame, inplace: bool = False):
        # synthesize a new column for obs_label_df by concatenating
        # the self.batch_column_names columns
        if not inplace:
            obs_label_df = obs_label_df.copy()
        obs_label_df[self.label_colname] = (
            obs_label_df[self.label_keys].astype(str).agg(self.labels_colsep.join, axis=1)
        )
        return obs_label_df

    def _add_categ_cov_col(self, obs_categ_cov_df: pd.DataFrame, inplace: bool = False):
        # synthesize a new column for obs_label_df by concatenating
        # the self.batch_column_names columns
        if not inplace:
            obs_categ_cov_df = obs_categ_cov_df.copy()
        obs_categ_cov_df[self._categorical_covariate_colname] = (
            obs_categ_cov_df[self._categorical_covariate_keys]
            .astype(str)
            .agg(self.categ_cov_colsep.join, axis=1)
        )
        return obs_categ_cov_df

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
            ).unsqueeze(1)
            if self.batch_column_names is not None
            else None,
            "labels": torch.from_numpy(
                self.label_encoder.transform(batch_obs[self.label_colname])
            ).unsqueeze(1)
            if self.label_keys is not None
            else torch.empty(0),
            "extra_categorical_covs": torch.cat(
                [
                    torch.from_numpy(self.categ_cov_encoder.transform(batch_obs[k])).unsqueeze(1)
                    for k in self._categorical_covariate_keys
                ],
                dim=1,
            )
            if self._categorical_covariate_keys is not None
            else None,
            "extra_continuous_covs": torch.cat(
                [
                    torch.from_numpy(batch_obs[k].values).float().unsqueeze(1)
                    for k in self._continuous_covariate_keys
                ],
                dim=1,
            )
            if self._continuous_covariate_keys is not None
            else None,
        }

    # scVI code expects these properties on the DataModule:

    @property
    def unlabeled_category(self) -> str:
        """String assigned to unlabeled cells."""
        if not hasattr(self, "_unlabeled_category"):
            raise AttributeError("`unlabeled_category` not set.")
        return self._unlabeled_category

    @unlabeled_category.setter
    def unlabeled_category(self, value: str | None):
        if not (value is None or isinstance(value, str)):
            raise ValueError("`unlabeled_category` must be a string or None.")
        self._unlabeled_category = value

    @property
    def split_seed(self) -> int:
        """Seed for data split."""
        if not hasattr(self, "_split_seed"):
            raise AttributeError("`split_seed` not set.")
        return self._split_seed

    @split_seed.setter
    def split_seed(self, value: int | None):
        if value is not None and not isinstance(value, int):
            raise ValueError("`split_seed` must be an integer.")
        self._split_seed = value or 0

    @property
    def train_size(self) -> float:
        """Fraction of data to use for training."""
        if not hasattr(self, "_train_size"):
            raise AttributeError("`train_size` not set.")
        return self._train_size

    @train_size.setter
    def train_size(self, value: float | None):
        if value is not None and not isinstance(value, float):
            raise ValueError("`train_size` must be a float.")
        elif value is not None and (value < 0.0 or value > 1.0):
            raise ValueError("`train_size` must be between 0.0 and 1.0.")
        self._train_size = value or 1.0

    @property
    def validation_size(self) -> float:
        """Fraction of data to use for validation."""
        if not hasattr(self, "_train_size"):
            raise AttributeError("`validation_size` not available.")
        return 1.0 - self.train_size

    @property
    def n_obs(self) -> int:
        return len(self.query.obs_joinids())

    @property
    def n_vars(self) -> int:
        return len(self.query.var_joinids())

    @property
    def n_batch(self) -> int:
        return len(self.batch_encoder.classes_)

    @property
    def n_labels(self) -> int:
        if self.label_keys is not None:
            return len(self.labels_mapping)
        else:
            return 0

    @property
    def labels_mapping(self) -> int:
        if self.label_keys is not None:
            combined = np.concatenate((self.label_encoder.classes_, [self.unlabeled_category]))
            unique_values, idx = np.unique(combined, return_index=True)
            unique_values = unique_values[np.argsort(idx)]
            return unique_values.astype(object)

    @property
    def extra_categorical_covs(self) -> dict:
        if self._categorical_covariate_keys is None:
            out = {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_categorical_covs": 0},
            }
        else:
            # TODO need to adjust this mapping
            mapping = dict(
                zip(
                    self._categorical_covariate_keys,
                    [self.categ_cov_encoder.classes_],
                    strict=False,
                )
            )
            out = {
                "data_registry": {"attr_key": "_scvi_extra_categorical_covs", "attr_name": "obsm"},
                "state_registry": {
                    "field_keys": self._categorical_covariate_keys,
                    "mapping": mapping,
                    "n_cats_per_key": [len(mapping[map]) for map in mapping.keys()],
                },
                "summary_stats": {
                    "n_extra_categorical_covs": len(self._categorical_covariate_keys)
                },
            }
        return out

    @property
    def extra_continuous_covs(self) -> dict:
        if self._continuous_covariate_keys is None:
            out = {
                "data_registry": {},
                "state_registry": {},
                "summary_stats": {"n_extra_continuous_covs": 0},
            }
        else:
            out = {
                "data_registry": {"attr_key": "_scvi_extra_continuous_covs", "attr_name": "obsm"},
                "state_registry": {
                    "columns": np.array(self._continuous_covariate_keys, dtype=object)
                },
                "summary_stats": {"n_extra_continuous_covs": len(self._continuous_covariate_keys)},
            }
        return out

    @property
    def registry(self) -> dict:
        features_names = list(
            self.query.var_joinids().tolist() if self.query is not None else range(self.n_vars)
        )
        return {
            "scvi_version": scvi.__version__,
            "model_name": self.model_name,
            "setup_args": {
                "layer": None,
                "batch_key": self.batch_colname,
                "labels_key": self.label_keys[0] if self.label_keys is not None else "label",
                "size_factor_key": None,
                "categorical_covariate_keys": self._categorical_covariate_keys,
                "continuous_covariate_keys": self._continuous_covariate_keys,
            },
            "field_registries": {
                "X": {
                    "data_registry": {"attr_name": "X", "attr_key": None},
                    "state_registry": {
                        "n_obs": self.n_obs,
                        "n_vars": self.n_vars,
                        "column_names": [str(i) for i in features_names],
                    },
                    "summary_stats": {"n_vars": self.n_vars, "n_cells": self.n_obs},
                },
                "batch": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_batch"},
                    "state_registry": {
                        "categorical_mapping": self.batch_labels,
                        "original_key": "batch",
                    },
                    "summary_stats": {"n_batch": self.n_batch},
                },
                "labels": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_labels"},
                    "state_registry": {
                        "categorical_mapping": self.labels_mapping,
                        "original_key": self.label_keys[0]
                        if self.label_keys is not None
                        else "label",
                        "unlabeled_category": self.unlabeled_category,
                    },
                    "summary_stats": {"n_labels": self.n_labels},
                },
                "size_factor": {"data_registry": {}, "state_registry": {}, "summary_stats": {}},
                "extra_categorical_covs": self.extra_categorical_covs,
                "extra_continuous_covs": self.extra_continuous_covs,
            },
            "setup_method_name": "setup_datamodule",
        }

    def inference_dataloader(self):
        """Dataloader for inference with `on_before_batch_transfer` applied."""
        dataloader = self.train_dataloader()
        return self._InferenceDataloader(dataloader, self.on_before_batch_transfer)

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
