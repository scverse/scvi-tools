from __future__ import annotations

import functools
from typing import TYPE_CHECKING

import numpy as np
import psutil
import torch
from cellxgene_census.experimental.ml import Encoder, ExperimentDataPipe, experiment_dataloader
from lightning.pytorch import LightningDataModule
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader

import scvi
from scvi.utils import dependencies

if TYPE_CHECKING:
    from typing import Any

    import numpy.typing as npt
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
class SCVIDataModule(LightningDataModule):
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
        self.batch_column_names = batch_column_names
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

        # dataset_kwargs = self.dataset_kwargs
        # if 'batch_keys' in dataset_kwargs.keys():
        #    dataset_kwargs.pop('batch_keys')
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


class BatchEncoder(Encoder):
    """An encoder that concatenates and encodes several obs columns."""

    def __init__(self, cols: list[str], name: str = "batch"):
        self.cols = cols
        from sklearn.preprocessing import LabelEncoder

        self._name = name
        self._encoder = LabelEncoder()

    def _join_cols(self, df: pd.DataFrame) -> pd.Series[str]:
        return functools.reduce(lambda a, b: a + b, [df[c].astype(str) for c in self.cols])

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs DataFrame into a DataFrame of encoded values."""
        arr = self._join_cols(df)
        return self._encoder.transform(arr)  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with obs."""
        arr = self._join_cols(obs)
        self._encoder.fit(arr.unique())

    @property
    def columns(self) -> list[str]:
        """Columns in `obs` that the encoder will be applied to."""
        return self.cols

    @property
    def name(self) -> str:
        """Name of the encoder."""
        return self._name

    @property
    def classes_(self) -> list[str]:
        """Classes of the encoder."""
        return self._encoder.classes_


class LabelEncoderNew(BatchEncoder):
    """An encoder that concatenates and encodes several obs columns as label +

    uses a string as missing observation label.
    """

    def __init__(self, cols: list[str], unlabeled_category: str = "Unknown", name: str = "label"):
        super().__init__(cols, name)
        self._unlabeled_category = unlabeled_category

    @property
    def unlabeled_category(self) -> str:
        """Name of the unlabeled_category."""
        return self._unlabeled_category

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs DataFrame into a DataFrame of encoded values + unlabeled category."""
        arr = self._join_cols(df)
        for unique_item in arr.unique():
            if unique_item not in self._encoder.classes_:
                arr = [self.unlabeled_category if x == unique_item else x for x in arr]
        return self._encoder.transform(arr)  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with obs + unlabeled category."""
        arr = self._join_cols(obs)
        self._encoder.fit(np.append(arr.unique(), self.unlabeled_category))


class CensusSCVIDataModule(LightningDataModule):
    """Lightning data module for training an scVI model using the ExperimentDataPipe.

    Parameters
    ----------
    *args
        Positional arguments passed to
        :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
    batch_keys
        List of obs column names concatenated to form the batch column.
    label_keys
        List of obs column names concatenated to form the label column.
    unlabeled_category
        Value used for unlabeled cells in `labels_key` used to set up CZI datamodule with scvi.
    train_size
        Fraction of data to use for training.
    split_seed
        Seed for data split.
    dataloader_kwargs
        Keyword arguments passed into
        :func:`~cellxgene_census.experimental.ml.pytorch.experiment_dataloader`.
    **kwargs
        Additional keyword arguments passed into
        :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`. Must not include
        ``obs_column_names``.
    """

    _TRAIN_KEY = "train"
    _VALIDATION_KEY = "validation"

    def __init__(
        self,
        *args,
        batch_keys: list[str] | None = None,
        label_keys: list[str] | None = None,
        unlabeled_category: str | None = "Unknown",
        train_size: float | None = None,
        split_seed: int | None = None,
        dataloader_kwargs: dict[str, any] | None = None,
        **kwargs,
    ):
        super().__init__()
        self.datapipe_args = args
        self.datapipe_kwargs = kwargs
        self.batch_keys = batch_keys
        self.label_keys = label_keys
        self.unlabeled_category = unlabeled_category
        self.train_size = train_size
        self.split_seed = split_seed
        self.dataloader_kwargs = dataloader_kwargs or {}

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
    def label_keys(self) -> list[str]:
        """List of obs column names concatenated to form the label column."""
        if not hasattr(self, "_label_keys"):
            raise AttributeError("`label_keys` not set.")
        return self._label_keys

    @label_keys.setter
    def label_keys(self, value: list[str] | None):
        if not (value is None or isinstance(value, list)):
            raise ValueError("`label_keys` must be a list of strings or None.")
        self._label_keys = value

    @property
    def batch_keys(self) -> list[str]:
        """List of obs column names concatenated to form the batch column."""
        if not hasattr(self, "_batch_keys"):
            raise AttributeError("`batch_keys` not set.")
        return self._batch_keys

    @batch_keys.setter
    def batch_keys(self, value: list[str] | None):
        if value is None or not isinstance(value, list):
            raise ValueError("`batch_keys` must be a list of strings.")
        self._batch_keys = value

    @property
    def obs_column_names(self) -> list[str]:
        """Passed to :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`."""
        if hasattr(self, "_obs_column_names"):
            return self._obs_column_names

        obs_column_names = []
        if self.batch_keys is not None:
            obs_column_names.extend(self.batch_keys)

        self._obs_column_names = obs_column_names
        return self._obs_column_names

    @property
    def obs_label_names(self) -> list[str]:
        """Passed to :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`."""
        if hasattr(self, "_obs_label_names"):
            return self._obs_label_names

        obs_label_names = []
        if self.label_keys is not None:
            obs_label_names.extend(self.label_keys)

        self._obs_label_names = obs_label_names
        return self._obs_label_names

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
    def weights(self) -> dict[str, float]:
        """Passed to :meth:`~cellxgene_census.experimental.ml.ExperimentDataPipe.random_split`."""
        if not hasattr(self, "_weights"):
            self._weights = {self._TRAIN_KEY: self.train_size}
            if self.validation_size > 0.0:
                self._weights[self._VALIDATION_KEY] = self.validation_size
        return self._weights

    @property
    def datapipe(self) -> ExperimentDataPipe:
        """Experiment data pipe."""
        if not hasattr(self, "_datapipe"):
            batch_encoder = BatchEncoder(self.obs_column_names)
            encoders_list = [batch_encoder]
            if self.label_keys is not None:
                label_encoder = LabelEncoderNew(self.obs_label_names, self.unlabeled_category)
                encoders_list.append(label_encoder)
            self._datapipe = ExperimentDataPipe(
                *self.datapipe_args,
                encoders=encoders_list,
                **self.datapipe_kwargs,
            )
        return self._datapipe

    def setup(self, stage: str | None = None):
        """Set up the train and validation data pipes."""
        datapipes = self.datapipe.random_split(weights=self.weights, seed=self.split_seed)
        self._train_datapipe = datapipes[0]
        if self.validation_size > 0.0:
            self._validation_datapipe = datapipes[1]
        else:
            self._validation_datapipe = None

    def train_dataloader(self):
        """Training data loader."""
        return experiment_dataloader(self._train_datapipe, **self.dataloader_kwargs)

    def val_dataloader(self):
        """Validation data loader."""
        if self._validation_datapipe is not None:
            return experiment_dataloader(self._validation_datapipe, **self.dataloader_kwargs)

    @property
    def n_obs(self) -> int:
        """Number of observations in the query.

        Necessary in scvi-tools to compute a heuristic of ``max_epochs``.
        """
        return self.datapipe.shape[0]

    @property
    def n_vars(self) -> int:
        """Number of features in the query.

        Necessary in scvi-tools to initialize the actual layers in the model.
        """
        return self.datapipe.shape[1]

    @property
    def n_batch(self) -> int:
        """Number of unique batches (after concatenation of ``batch_keys``).

        Necessary in scvi-tools so that the model knows how to one-hot encode batches.
        """
        return self.get_n_classes("batch")

    @property
    def n_label(self) -> int:
        """Number of unique labels (after concatenation of ``label_keys``).

        Necessary in scvi-tools so that the model knows how to one-hot encode labels.
        """
        return self.get_n_classes("label")

    def get_n_classes(self, key: str) -> int:
        """Return the number of classes for a given obs column."""
        return len(self.datapipe.obs_encoders[key].classes_)

    def on_before_batch_transfer(
        self,
        batch: tuple[torch.Tensor, torch.Tensor],
        dataloader_idx: int,
    ) -> dict[str, torch.Tensor | None]:
        """Format the datapipe output with registry keys for scvi-tools."""
        X, obs = batch

        X_KEY: str = "X"
        BATCH_KEY: str = "batch"
        LABELS_KEY: str = "labels"

        return {
            X_KEY: X,
            BATCH_KEY: obs,
            LABELS_KEY: torch.empty(0),  # None,
        }
