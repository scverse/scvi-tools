from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from sklearn.cluster import KMeans

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type, get_anndata_attribute
from scvi.data.fields import CategoricalJointObsField, CategoricalObsField, LayerField
from scvi.model.base import (
    BaseMinifiedModeModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module import VAEC
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseMinifiedModeModelClass):
    """Conditional version of single-cell Variational Inference, used for multi-resolution deconvolution of spatial transcriptomics data :cite:p:`Lopez22`.

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
        Whether to reweight observations by their inverse proportion (useful for lowly abundant cell types)
    dropout_rate
        Dropout rate for neural networks.
    **module_kwargs
        Keyword args for :class:`~scvi.modules.VAEC`

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
    _LATENT_QZM = "_condscvi_latent_qzm"
    _LATENT_QZV = "_condscvi_latent_qzv"
    _OBSERVED_LIB_SIZE = "_condscvi_observed_lib_size"

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
        if 'fine_labels' in self.adata_manager.data_registry:
            fine_labels = get_anndata_attribute(
                adata,
                self.adata_manager.data_registry.labels.attr_name,
                '_scvi_fine_labels',
            )
            coarse_labels = get_anndata_attribute(
                adata,
                self.adata_manager.data_registry.labels.attr_name,
                '_scvi_labels'
            )

            df_ct = pd.DataFrame({
                'fine_labels_key': fine_labels.ravel(),
                'coarse_labels_key': coarse_labels.ravel()}).drop_duplicates()
            print('YYYY', df_ct)
            fine_labels_mapping = self.adata_manager.get_state_registry(
                'fine_labels'
            ).categorical_mapping
            coarse_labels_mapping = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
            ).categorical_mapping

            df_ct['fine_labels'] = fine_labels_mapping[df_ct['fine_labels_key']]
            df_ct['coarse_labels'] = coarse_labels_mapping[df_ct['coarse_labels_key']]

            self.df_ct_name_dict = {}
            self.df_ct_id_dict = {}
            for i, row in df_ct.iterrows():
                count = len(df_ct.loc[:i][df_ct['coarse_labels'] == row['coarse_labels']]) - 1
                self.df_ct_name_dict[row['fine_labels']] = (row['coarse_labels'], row['coarse_labels_key'], count)
                self.df_ct_id_dict[row['fine_labels_key']] = (row['coarse_labels'], row['coarse_labels_key'], count)
        else:
            self.df_ct_name_dict = None
            self.df_ct_id_dict = None

        self.module = self._module_cls(
            n_input=n_vars,
            n_batch=n_batch,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            df_ct_id_dict=self.df_ct_id_dict,
            **module_kwargs,
        )
        self._model_summary_string = (
            "Conditional SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: {}, weight_obs: {}"
        ).format(n_hidden, n_latent, n_layers, dropout_rate, weight_obs)
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def get_vamp_prior(self, adata: AnnData | None = None, p: int = 10, scales_prior_n_samples: int | None = None, default_cat: list | None = None) -> np.ndarray:
        r"""Return an empirical prior over the cell-type specific latent space (vamp prior) that may be used for deconvolution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        p
            number of clusters in kmeans clustering for cell-type sub-clustering for empirical prior
        scales_prior_n_samples
            return scales of negative binomial distribution for calculates prior means and variances using n_samples.
        default_cat
            default value for categorical covariates

        Returns
        -------
        mean_vprior: np.ndarray
            (n_labels, p, D) array
        var_vprior
            (n_labels, p, D) array
        weights_vprior
            (n_labels, p) array
        scales_vprior
            (n_labels, p, G) array
        """
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
                "var_vprior": torch.exp(self.module.prior_log_scales)**2 + 1e-4,
                "weights_vprior": torch.nn.functional.softmax(self.module.prior_logits, dim=-1)
            }
        else:

            # Extracting latent representation of adata including variances.
            mean_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
            var_vprior = np.ones((self.summary_stats.n_labels, p, self.module.n_latent))
            mp_vprior = np.zeros((self.summary_stats.n_labels, p))

            labels_state_registry = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.LABELS_KEY
            )
            key = labels_state_registry.original_key
            mapping = labels_state_registry.categorical_mapping

            mean_cat, var_cat = self.get_latent_representation(adata, return_dist=True)

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
                        Given cell type specific clustering contains more clusters than vamp_prior_p.
                        Increase value of vamp_prior_p to largest number of cell type specific clusters."""

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
                "weights_vprior": mp_vprior
            }

        if scales_prior_n_samples is not None:
            scales_vprior = np.zeros((self.summary_stats.n_labels, p, self.summary_stats.n_vars))
            cat_covs = [
                torch.full([scales_prior_n_samples, 1], float(np.where(value==default_cat[ind])[0]) if default_cat else 0, device=self.module.device)
                           for ind, value in enumerate(self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).mappings.values())]
            for ct in range(self.summary_stats["n_labels"]):
                for cluster in range(p):
                    sampled_z = torch.distributions.Normal(
                        results['mean_vprior'][ct, cluster, :], torch.sqrt(results['var_vprior'][ct, cluster, :])
                    ).sample([scales_prior_n_samples,]).to(self.module.device)
                    h = self.module.decoder(sampled_z, torch.full([scales_prior_n_samples, 1], ct, device=self.module.device), *cat_covs)
                    scales_vprior[ct, cluster, :] = self.module.px_decoder(h).mean(0).cpu()
            results["scales_vprior"] = scales_vprior

        return results

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
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set are split in the
            sequential order of the data according to `validation_size` and `train_size` percentages.
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

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        fine_labels_key: str | None = None,
        layer: str | None = None,
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
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        if fine_labels_key is not None:
            anndata_fields.append(
                CategoricalObsField('fine_labels', fine_labels_key
                )
            )

        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(cls, adata_minify_type)
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
