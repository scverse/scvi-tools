import logging
import warnings
from typing import Optional, Union

import numpy as np
import torch
from anndata import AnnData
import scanpy
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
    Conditional version of single-cell Variational Inference, used for multi-resolution deconvolution of spatial transcriptomics data [Lopez21]_.

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
    dropout_rate
        Dropout rate for the encoder neural networks.
    weight_obs
        Whether to reweight observations by their inverse proportion (useful for lowly abundant cell types)
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

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        weight_obs: bool = False,
        **module_kwargs,
    ):
        super(CondSCVI, self).__init__(adata)

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

        self.module = VAEC(
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

    @torch.no_grad()
    def get_vamp_prior(
        self,
        adata: Optional[AnnData] = None,
        p: int = 50,
        resolution: float = 10.0
    ) -> np.ndarray:
        r"""
        Return an empirical prior over the cell-type specific latent space (vamp prior) that may be used for deconvolution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        p
            number of components in the mixture model underlying the empirical prior
        resolution
            resolution of overclustering in leiden used for metapoints in empirical prior

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

        mean_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
        var_vprior = np.zeros((self.summary_stats.n_labels, p, self.module.n_latent))
        labels_state_registry = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        key = labels_state_registry.original_key
        mapping = labels_state_registry.categorical_mapping

        scdl = self._make_data_loader(
            adata=adata, batch_size=p
        )

        mean = []
        var = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            y = tensors[REGISTRY_KEYS.LABELS_KEY]
            out = self.module.inference(x, y)
            mean_, var_ = out["qz_m"], out["qz_v"]
            mean += [mean_.cpu()]
            var += [var_.cpu()]

        mean_cat, var_cat = np.array(torch.cat(mean)), np.array(torch.cat(var))
        adata.obsm['X_CondSCVI'] = mean_cat

        for ct in range(self.summary_stats["n_labels"]):
            local_indices = np.where(adata.obs[key] == mapping[ct])[0]
            sub_adata = adata[local_indices, :].copy()
            scanpy.pp.neighbors(sub_adata, use_rep="X_CondSCVI")
            if "overclustering_vamp" not in sub_adata.obs.columns:
                scanpy.tl.leiden(sub_adata, resolution=resolution, key_added="overclustering_vamp")

            var_cluster = np.zeros([len(np.unique(sub_adata.obs.overclustering_vamp)), self.module.n_latent])
            mean_cluster = np.zeros([len(np.unique(sub_adata.obs.overclustering_vamp)), self.module.n_latent])
            

            for j in np.unique(sub_adata.obs.overclustering_vamp):
                indices_curr = local_indices[np.where(sub_adata.obs['overclustering_vamp'] == j)[0]]
                var_cluster[int(j),:] = np.mean(var_cat[indices_curr], axis=0) + \
                                        np.mean(mean_cat[indices_curr]**2, axis=0) - \
                                        (np.mean(mean_cat[indices_curr], axis=0))**2
                mean_cluster[int(j),:] = np.mean(mean_cat[indices_curr], axis=0)

            unique, counts = np.unique(sub_adata.obs.overclustering_vamp, return_counts=True)
            selection = np.random.choice(np.arange(len(unique)), size=p, p=counts/sum(counts))
            mean_vprior[ct, :, :] = mean_cluster[selection, :]
            var_vprior[ct, :, :] = var_cluster[selection, :]
        
        return mean_vprior, var_vprior

    def train(
        self,
        max_epochs: int = 400,
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
