import logging
import warnings
from typing import Optional, Union

import numpy as np
import torch
from anndata import AnnData
from sklearn.cluster import KMeans

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.model.base import (
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module import VAEC
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """
    Conditional version of single-cell Variational Inference, used for multi-resolution deconvolution of spatial transcriptomics data :cite:p:`Lopez21`.

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

        self.module = self._module_cls(
            n_input=n_vars,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **module_kwargs,
        )
        self._model_summary_string = (
            "Conditional SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: {}, weight_obs: {}"
        ).format(n_hidden, n_latent, n_layers, dropout_rate, weight_obs)
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def get_vamp_prior(
        self, adata: Optional[AnnData] = None, p: int = 10
    ) -> np.ndarray:
        r"""
        Return an empirical prior over the cell-type specific latent space (vamp prior) that may be used for deconvolution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        p
            number of clusters in kmeans clustering for cell-type sub-clustering for empirical prior

        Returns
        -------
        mean_vprior: np.ndarray
            (n_labels, p, D) array
        var_vprior
            (n_labels, p, D) array
        """
        if self.is_trained_ is False:
            warnings.warn(
                "Trying to query inferred values from an untrained model. Please train the model first."
            )

        adata = self._validate_anndata(adata)

        # Extracting latent representation of adata including variances.
        mean_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
        var_vprior = np.ones((self.summary_stats.n_labels, p, self.module.n_latent))
        mp_vprior = np.zeros((self.summary_stats.n_labels, p))

        labels_state_registry = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        key = labels_state_registry.original_key
        mapping = labels_state_registry.categorical_mapping

        scdl = self._make_data_loader(adata=adata, batch_size=p)

        mean = []
        var = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            y = tensors[REGISTRY_KEYS.LABELS_KEY]
            out = self.module.inference(x, y)
            mean_, var_ = out["qz"].loc, (out["qz"].scale ** 2)
            mean += [mean_.cpu()]
            var += [var_.cpu()]

        mean_cat, var_cat = torch.cat(mean).numpy(), torch.cat(var).numpy()

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
                indices_curr = local_indices[
                    np.where(overclustering_vamp == cluster)[0]
                ]
                var_cluster[index, :] = np.mean(var_cat[indices_curr], axis=0) + np.var(
                    mean_cat[indices_curr], axis=0
                )
                mean_cluster[index, :] = np.mean(mean_cat[indices_curr], axis=0)

            slicing = slice(n_labels_overclustering)
            mean_vprior[ct, slicing, :] = mean_cluster
            var_vprior[ct, slicing, :] = var_cluster
            mp_vprior[ct, slicing] = counts / sum(counts)

        return mean_vprior, var_vprior, mp_vprior

    def train(
        self,
        max_epochs: int = 300,
        lr: float = 0.001,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using MAP inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
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
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_labels_key)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
