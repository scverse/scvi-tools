from __future__ import annotations

import importlib
import logging
from typing import TYPE_CHECKING

import anndata
import numpy as np
import pandas as pd
import torch

from scvi import REGISTRY_KEYS
from scvi.data._utils import get_anndata_attribute
from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter
from scvi.model._utils import get_max_epochs_heuristic, use_distributed_sampler
from scvi.train import SemiSupervisedTrainingPlan, TrainingPlan, TrainRunner
from scvi.train._callbacks import SubSampleLabels
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from collections.abc import Sequence

    from lightning import LightningDataModule

    from scvi._types import AnnOrMuData

logger = logging.getLogger(__name__)


class UnsupervisedTrainingMixin:
    """General purpose unsupervised train method."""

    _data_splitter_cls = DataSplitter
    _training_plan_cls = TrainingPlan
    _train_runner_cls = TrainRunner

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        load_sparse_tensor: bool = False,
        batch_size: int = 128,
        early_stopping: bool = False,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        datamodule: LightningDataModule | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            The maximum number of epochs to train the model. The actual number of epochs may be
            less if early stopping is enabled. If ``None``, defaults to a heuristic based on
            :func:`~scvi.model.get_max_epochs_heuristic`. Must be passed in if ``datamodule`` is
            passed in, and it does not have an ``n_obs`` attribute.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range ``[0.0, 1.0]``. Passed into
            :class:`~scvi.dataloaders.DataSplitter`. Not used if ``datamodule`` is passed in.
        validation_size
            Size of the test set. If ``None``, defaults to ``1 - train_size``. If
            ``train_size + validation_size < 1``, the remaining cells belong to a test set. Passed
            into :class:`~scvi.dataloaders.DataSplitter`. Not used if ``datamodule`` is passed in.
        shuffle_set_split
            Whether to shuffle indices before splitting. If ``False``, the val, train, and test set
            are split in the sequential order of the data according to ``validation_size`` and
            ``train_size`` percentages. Passed into :class:`~scvi.dataloaders.DataSplitter`. Not
            used if ``datamodule`` is passed in.
        load_sparse_tensor
            ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
            :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to
            GPUs, depending on the sparsity of the data. Passed into
            :class:`~scvi.dataloaders.DataSplitter`. Not used if ``datamodule`` is passed in.
        batch_size
            Minibatch size to use during training. Passed into
            :class:`~scvi.dataloaders.DataSplitter`. Not used if ``datamodule`` is passed in.
        early_stopping
            Perform early stopping. Additional arguments can be passed in through ``**kwargs``.
            See :class:`~scvi.train.Trainer` for further options.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
            Values in this argument can be overwritten by arguments directly passed into this
            method, when appropriate. Not used if ``datamodule`` is passed in.
        plan_kwargs
            Additional keyword arguments passed into :class:`~scvi.train.TrainingPlan`. Values in
            this argument can be overwritten by arguments directly passed into this method, when
            appropriate.
        datamodule
            ``EXPERIMENTAL`` A :class:`~lightning.pytorch.core.LightningDataModule` instance to use
            for training in place of the default :class:`~scvi.dataloaders.DataSplitter`. Can only
            be passed in if the model was not initialized with :class:`~anndata.AnnData`.
        **kwargs
           Additional keyword arguments passed into :class:`~scvi.train.Trainer`.
        """
        if datamodule is not None and not self._module_init_on_train:
            raise ValueError(
                "Cannot pass in `datamodule` if the model was initialized with `adata`."
            )
        elif datamodule is None and self._module_init_on_train:
            raise ValueError(
                "If the model was not initialized with `adata`, a `datamodule` must be passed in."
            )

        if max_epochs is None:
            if datamodule is None:
                max_epochs = get_max_epochs_heuristic(self.adata.n_obs)
            elif hasattr(datamodule, "n_obs"):
                max_epochs = get_max_epochs_heuristic(datamodule.n_obs)
            else:
                raise ValueError(
                    "If `datamodule` does not have `n_obs` attribute, `max_epochs` must be "
                    "passed in."
                )

        if datamodule is None:
            datasplitter_kwargs = datasplitter_kwargs or {}
            datamodule = self._data_splitter_cls(
                self.adata_manager,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
                distributed_sampler=use_distributed_sampler(trainer_kwargs.get("strategy", None)),
                load_sparse_tensor=load_sparse_tensor,
                **datasplitter_kwargs,
            )
        elif self.module is None:
            self.module = self._module_cls(
                datamodule.n_vars,
                n_batch=datamodule.n_batch,
                n_labels=getattr(datamodule, "n_labels", 1),
                n_continuous_cov=getattr(datamodule, "n_continuous_cov", 0),
                n_cats_per_cov=getattr(datamodule, "n_cats_per_cov", None),
                **self._module_kwargs,
            )

        plan_kwargs = plan_kwargs or {}
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=datamodule,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            **trainer_kwargs,
        )
        return runner()


class SemisupervisedTrainingMixin:
    _training_plan_cls = SemiSupervisedTrainingPlan

    def _set_indices_and_labels(self):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        self.original_label_key = labels_state_registry.original_key
        self.unlabeled_category_ = labels_state_registry.unlabeled_category

        labels = get_anndata_attribute(
            self.adata,
            self.adata_manager.data_registry.labels.attr_name,
            self.original_label_key,
            mod_key=getattr(self.adata_manager.data_registry.labels, "mod_key", None),
        ).ravel()
        self._label_mapping = labels_state_registry.categorical_mapping

        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category_).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category_).ravel()
        self._code_to_label = dict(enumerate(self._label_mapping))

    def predict(
        self,
        adata: AnnOrMuData | None = None,
        indices: Sequence[int] | None = None,
        soft: bool = False,
        batch_size: int | None = None,
        use_posterior_mean: bool = True,
        ig_interpretability: bool = False,
        ig_args: dict | None = None,
    ) -> (np.ndarray | pd.DataFrame, None | np.ndarray):
        """Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData or MuData object that has been registered via corresponding setup
            method in model class.
        indices
            Return probabilities for each class label.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        use_posterior_mean
            If ``True``, uses the mean of the posterior distribution to predict celltype
            labels. Otherwise, uses a sample from the posterior distribution - this
            means that the predictions will be stochastic.
        ig_interpretability
            If True, run the integrated circuits interpretability per sample and returns a score
            matrix, in which for each sample we score each gene for its contribution to the
            sample prediction
        ig_args
            Keyword args for IntegratedGradients
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        attributions = None
        if ig_interpretability:
            missing_modules = []
            try:
                importlib.import_module("captum")
            except ImportError:
                missing_modules.append("captum")
            if len(missing_modules) > 0:
                raise ModuleNotFoundError("Please install captum to use this functionality.")
            from captum.attr import IntegratedGradients

            ig = IntegratedGradients(self.module.classify)
            attributions = []

        # in case of no indices to predict return empty values
        if len(indices) == 0:
            pred = []
            if ig_interpretability:
                return pred, attributions
            else:
                return pred

        scdl = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )

        y_pred = []
        for _, tensors in enumerate(scdl):
            inference_inputs = self.module._get_inference_input(tensors)  # (n_obs, n_vars)
            data_inputs = {
                key: inference_inputs[key]
                for key in inference_inputs.keys()
                if key not in ["batch_index", "cont_covs", "cat_covs"]
            }

            batch = tensors[REGISTRY_KEYS.BATCH_KEY]

            cont_key = REGISTRY_KEYS.CONT_COVS_KEY
            cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

            cat_key = REGISTRY_KEYS.CAT_COVS_KEY
            cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

            pred = self.module.classify(
                **data_inputs,
                batch_index=batch,
                cat_covs=cat_covs,
                cont_covs=cont_covs,
                use_posterior_mean=use_posterior_mean,
            )
            if self.module.classifier.logits:
                pred = torch.nn.functional.softmax(pred, dim=-1)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

            if ig_interpretability:
                # we need the hard prediction if was not done yet
                hard_pred = pred.argmax(dim=1) if soft else pred
                ig_args = ig_args or {}
                attribution = ig.attribute(
                    tuple(data_inputs.values()), target=hard_pred, **ig_args
                )
                attributions.append(attribution[0])

        if ig_interpretability:
            if attributions is not None and len(attributions) > 0:
                attributions = torch.cat(attributions, dim=0).detach().numpy()
                attributions = self.get_ranked_genes(adata, attributions)

        if len(y_pred) > 0:
            y_pred = torch.cat(y_pred).numpy()
            if not soft:
                predictions = [self._code_to_label[p] for p in y_pred]
                if ig_interpretability:
                    return np.array(predictions), attributions
                else:
                    return np.array(predictions)
            else:
                n_labels = len(pred[0])
                pred = pd.DataFrame(
                    y_pred,
                    columns=self._label_mapping[:n_labels],
                    index=adata.obs_names[indices],
                )
                if ig_interpretability:
                    return pred, attributions
                else:
                    return pred

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        n_samples_per_label: float | None = None,
        check_val_every_n_epoch: int | None = None,
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train,
            and test set are split in the sequential order of the data according to
            `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        %(param_accelerator)s
        %(param_devices)s
        datasplitter_kwargs
            Additional keyword arguments passed into
            :class:`~scvi.dataloaders.SemiSupervisedDataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.SemiSupervisedTrainingPlan`. Keyword
            arguments passed to `train()` will overwrite values present in `plan_kwargs`,
            when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

            if self.was_pretrained:
                max_epochs = int(np.min([10, np.max([2, round(max_epochs / 3.0)])]))

        logger.info(f"Training for {max_epochs} epochs.")

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs
        datasplitter_kwargs = datasplitter_kwargs or {}

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = [SubSampleLabels()] if len(self._labeled_indices) != 0 else []

        data_splitter = SemiSupervisedDataSplitter(
            adata_manager=self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            n_samples_per_label=n_samples_per_label,
            distributed_sampler=use_distributed_sampler(trainer_kwargs.get("strategy", None)),
            batch_size=batch_size,
            **datasplitter_kwargs,
        )
        training_plan = self._training_plan_cls(self.module, self.n_labels, **plan_kwargs)

        if "callbacks" in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] + [sampler_callback]
        else:
            trainer_kwargs["callbacks"] = sampler_callback

        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **trainer_kwargs,
        )
        return runner()

    def get_ranked_genes(
        self, adata: AnnOrMuData | None = None, attrs: np.ndarray | None = None
    ) -> pd.DataFrame:
        """Get the ranked gene list based on highest attributions.

        Parameters
        ----------
        attr: numpy.ndarray
            Attributions matrix.

        Returns
        -------
        pandas.DataFrame
            A pandas dataframe containing the ranked attributions for each gene

        Examples
        --------
        >>> attrs_df = interpreter.get_ranked_genes(attrs)
        """
        if attrs is None:
            Warning("Missing Attributions matrix")
            return

        adata = self._validate_anndata(adata)

        mean_attrs = attrs.mean(axis=0)
        idx = mean_attrs.argsort()[::-1]
        df = {
            "gene": np.array(adata.var_names)[idx],
            "gene_idx": idx,
            "attribution_mean": mean_attrs[idx],
            "attribution_std": attrs.std(axis=0)[idx],
            "cells": attrs.shape[0],
        }
        return pd.DataFrame(df)

    def shap_adata_predict(
        self,
        X,
    ):
        adata = self._validate_anndata()

        # we need to adjust adata to the shap random selection ..
        if len(X) > len(adata):
            # Repeat the data to expand to a larger size
            n_repeats = len(X) / len(adata)  # how many times you want to repeat the data
            adata_to_pred = adata[adata.obs.index.repeat(n_repeats), :]
            if len(X) > len(adata_to_pred):
                adata_to_pred = anndata.concat(
                    [adata_to_pred, adata[0 : (len(X) - len(adata_to_pred))]]
                )
        else:
            adata_to_pred = adata[0 : len(X)]
        adata_to_pred.X = X

        return self.predict(adata_to_pred, soft=True)

    def shap_predict(self, adata: AnnOrMuData | None = None, max_size: int = 100):
        missing_modules = []
        try:
            importlib.import_module("shap")
        except ImportError:
            missing_modules.append("shap")
        if len(missing_modules) > 0:
            raise ModuleNotFoundError("Please install shap to use this functionality.")
        import shap

        adata_orig = self._validate_anndata()
        adata = self._validate_anndata(adata)

        if type(adata_orig.X).__name__ == "csr_matrix":
            feature_matrix_background = pd.DataFrame.sparse.from_spmatrix(
                adata_orig.X, columns=adata_orig.var_names
            )
        else:
            feature_matrix_background = pd.DataFrame(adata_orig.X, columns=adata_orig.var_names)
        if type(adata.X).__name__ == "csr_matrix":
            feature_matrix = pd.DataFrame.sparse.from_spmatrix(
                adata.X, columns=adata_orig.var_names
            )
        else:
            feature_matrix = pd.DataFrame(adata.X, columns=adata_orig.var_names)
        feature_matrix_background = shap.sample(feature_matrix_background, max_size)
        feature_matrix = shap.sample(feature_matrix, max_size)
        explainer = shap.KernelExplainer(self.shap_adata_predict, feature_matrix_background)
        shap_values = explainer.shap_values(feature_matrix)
        return shap_values
