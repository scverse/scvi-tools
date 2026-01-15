from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalObsField,
)
from scvi.model.base import (
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module import VAEC
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """Conditional version of single-cell Variational Inference.

    Used for multi-resolution deconvolution of spatial transcriptomics data :cite:p:`Lopez22`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.CondSCVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    weight_obs
        Whether to reweight observations by their inverse proportion (useful for lowly abundant
        cell types)
    dropout_rate
        Dropout rate for neural networks.
    **module_kwargs
        Keyword args for :class:`~scvi.module.VAEC`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.CondSCVI.setup_anndata(adata, "labels")
    >>> vae = scvi.model.CondSCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_CondSCVI"] = vae.get_latent_representation()

    Notes
    -----
    See further usage examples in the following tutorial:

    1. :doc:`/tutorials/notebooks/spatial/DestVI_tutorial`
    """

    _module_cls = VAEC

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        weight_obs: bool = False,
        dropout_rate: float = 0.05,
        **module_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch
        n_labels = self.summary_stats.n_labels
        n_vars = self.summary_stats.n_vars
        if "n_fine_labels" in self.summary_stats:
            self.n_fine_labels = self.summary_stats.n_fine_labels
        else:
            self.n_fine_labels = None
        self._set_indices_and_labels(adata)
        if weight_obs:
            ct_counts = np.unique(
                self.get_from_registry(adata, REGISTRY_KEYS.LABELS_KEY),
                return_counts=True,
            )[1]
            ct_prop = ct_counts / np.sum(ct_counts)
            ct_prop[ct_prop < 0.05] = 0.05
            ct_prop = ct_prop / np.sum(ct_prop)
            ct_weight = 1.0 / ct_prop
            module_kwargs.update({"ct_weight": ct_weight})

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=getattr(self.summary_stats, "n_batch", 0),
            n_labels=self.summary_stats.n_labels,
            n_fine_labels=self.n_fine_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **module_kwargs,
        )
        self._model_summary_string = (
            f"Conditional SCVI Model with the following params: \nn_hidden: {n_hidden}, "
            f"n_latent: {n_latent}, n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"weight_obs: {weight_obs}"
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def get_vamp_prior(self, adata: AnnData | None = None, p: int = 10) -> np.ndarray:
        r"""Return an empirical prior over the cell-type specific latent space (vamp prior).

        May be used for deconvolution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        p
            number of clusters in kmeans clustering for cell-type sub-clustering for empirical
            prior

        Returns
        -------
        mean_vprior: np.ndarray
            (n_labels, p, D) array
        var_vprior
            (n_labels, p, D) array
        weights_vprior
            (n_labels, p) array
        """
        from sklearn.cluster import KMeans

        if self.is_trained_ is False:
            warnings.warn(
                "Trying to query inferred values from an untrained model. Please train "
                "the model first.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        adata = self._validate_anndata(adata)

        if self.module.prior == "mog":
            results = {
                "mean_vprior": self.module.prior_means,
                "var_vprior": torch.exp(self.module.prior_log_std) ** 2,
                "weights_vprior": torch.nn.functional.softmax(self.module.prior_logits, dim=-1),
            }
        else:
            # Extracting latent representation of adata including variances.
            mean_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
            var_vprior = np.ones((self.summary_stats.n_labels, p, self.module.n_latent))
            mp_vprior = np.zeros((self.summary_stats.n_labels, p))

            labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
            key = labels_state_registry.original_key
            mapping = labels_state_registry.categorical_mapping

            mean_cat, var_cat = self.get_latent_representation(
                adata, return_dist=True, batch_size=p, give_mean=False
            )

            for ct in range(self.summary_stats["n_labels"]):
                local_indices = np.where(adata.obs[key] == mapping[ct])[0]
                n_local_indices = len(local_indices)
                if "overclustering_vamp" not in adata.obs.columns:
                    if p < n_local_indices and p > 0:
                        overclustering_vamp = KMeans(n_clusters=p, n_init=30).fit_predict(
                            mean_cat[local_indices]
                        )
                    else:
                        # Every cell is its own cluster
                        overclustering_vamp = np.arange(n_local_indices)
                else:
                    overclustering_vamp = adata[local_indices, :].obs["overclustering_vamp"]

                keys, counts = np.unique(overclustering_vamp, return_counts=True)

                n_labels_overclustering = len(keys)
                if n_labels_overclustering > p:
                    error_mess = """
                        Given cell type specific clustering contains more clusters than
                        vamp_prior_p. Increase value of vamp_prior_p to largest number of cell type
                         specific clusters."""
                    raise ValueError(error_mess)

                var_cluster = np.ones(
                    [
                        n_labels_overclustering,
                        self.module.n_latent,
                    ]
                )
                mean_cluster = np.zeros_like(var_cluster)

                for index, cluster in enumerate(keys):
                    indices_curr = local_indices[np.where(overclustering_vamp == cluster)[0]]
                    var_cluster[index, :] = np.mean(var_cat[indices_curr], axis=0) + np.var(
                        mean_cat[indices_curr], axis=0
                    )
                    mean_cluster[index, :] = np.mean(mean_cat[indices_curr], axis=0)

                slicing = slice(n_labels_overclustering)
                mean_vprior[ct, slicing, :] = mean_cluster
                var_vprior[ct, slicing, :] = var_cluster
                mp_vprior[ct, slicing] = counts / sum(counts)
            results = {
                "mean_vprior": mean_vprior,
                "var_vprior": var_vprior,
                "weights_vprior": mp_vprior,
            }

        return results

    @torch.inference_mode()
    def predict(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        soft: bool = False,
        batch_size: int | None = None,
        use_posterior_mean: bool = True,
    ) -> np.ndarray | pd.DataFrame:
        """Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
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
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )
        y_pred = []
        for _, tensors in enumerate(scdl):
            inference_input = self.module._get_inference_input(tensors)
            qz = self.module.inference(**inference_input)["qz"]
            if use_posterior_mean:
                z = qz.loc
            else:
                z = qz.sample()
            pred = self.module.classify(
                z,
                label_index=inference_input["y"],
            )
            if self.module.classifier.logits:
                pred = torch.nn.functional.softmax(pred, dim=-1)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

        y_pred = torch.cat(y_pred).numpy()
        if not soft:
            predictions = []
            for p in y_pred:
                predictions.append(self._code_to_fine_label[p])

            return np.array(predictions)
        else:
            n_labels = len(pred[0])
            pred = pd.DataFrame(
                y_pred,
                columns=self._fine_label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
            return pred

    @torch.inference_mode()
    def confusion_coarse_celltypes(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
    ) -> np.ndarray | pd.DataFrame:
        """Return likelihood ratios of switching coarse cell-types

        to inform whether resolution is to granular.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
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
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )
        # Iterate once over the data and computes the reconstruction error
        keys = list(self._label_mapping) + ["original"]
        log_lkl = {key: [] for key in keys}
        for tensors in scdl:
            loss_kwargs = {"kl_weight": 1}
            _, _, losses = self.module(tensors, loss_kwargs=loss_kwargs)
            log_lkl["original"] += [losses.reconstruction_loss]
            for i in range(self.module.n_labels):
                tensors_ = tensors
                tensors_["y"] = torch.full_like(tensors["y"], i)
                _, _, losses = self.module(tensors_, loss_kwargs=loss_kwargs)
                log_lkl[keys[i]] += [losses.reconstruction_loss]
        for key in keys:
            log_lkl[key] = torch.stack(log_lkl[key]).detach().numpy()

        return log_lkl

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 300,
        lr: float = 0.001,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 1,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Trains the model using MAP inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            datasplitter_kwargs=datasplitter_kwargs,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    def _set_indices_and_labels(self, adata: AnnData):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        self.original_label_key = labels_state_registry.original_key
        self._label_mapping = labels_state_registry.categorical_mapping
        self._code_to_label = dict(enumerate(self._label_mapping))
        if self.n_fine_labels is not None:
            fine_labels_state_registry = self.adata_manager.get_state_registry("fine_labels")
            self.original_fine_label_key = fine_labels_state_registry.original_key
            self._fine_label_mapping = fine_labels_state_registry.categorical_mapping
            self._code_to_fine_label = dict(enumerate(self._fine_label_mapping))

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        fine_labels_key: str | None = None,
        layer: str | None = None,
        unlabeled_category: str = "unlabeled",
        size_factor_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        %(param_labels_key)s
        fine_labels_key
            Key in `adata.obs` where fine-grained labels are stored.
        %(param_layer)s
        %(param_unlabeled_category)s
        %(param_size_factor_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
        ]
        if fine_labels_key is not None:
            anndata_fields.append(
                LabelsWithUnlabeledObsField("fine_labels", fine_labels_key, unlabeled_category),
            )

        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
