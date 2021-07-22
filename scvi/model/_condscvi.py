import logging
import warnings
from typing import Optional, Union

import numpy as np
import torch
from anndata import AnnData

from scvi import _CONSTANTS
from scvi.model.base import (
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module import VAEC

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """
    Conditional version of single-cell Variational Inference, used for multi-resolution deconvolution of spatial transcriptomics data [Lopez21]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
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
    >>> scvi.data.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.external.CondSCVI(adata)
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

        n_labels = self.summary_stats["n_labels"]
        n_vars = self.summary_stats["n_vars"]
        if weight_obs:
            ct_counts = adata.obs["_scvi_labels"].value_counts()[range(n_labels)].values
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

        Returns
        -------
        mean_vprior: np.ndarray
            (n_labels, p, D) array
        var_vprior
            (n_labels, p, 3) array
        """
        if self.is_trained_ is False:
            warnings.warn(
                "Trying to query inferred values from an untrained model. Please train the model first."
            )

        adata = self._validate_anndata(adata)

        mean_vprior = np.zeros(
            (self.summary_stats["n_labels"], p, self.module.n_latent)
        )
        var_vprior = np.zeros((self.summary_stats["n_labels"], p, self.module.n_latent))
        key = self.scvi_setup_dict_["categorical_mappings"]["_scvi_labels"][
            "original_key"
        ]
        mapping = self.scvi_setup_dict_["categorical_mappings"]["_scvi_labels"][
            "mapping"
        ]
        for ct in range(self.summary_stats["n_labels"]):
            # pick p cells
            local_indices = np.random.choice(
                np.where(adata.obs[key] == mapping[ct])[0], p
            )
            # get mean and variance from posterior
            scdl = self._make_data_loader(
                adata=adata, indices=local_indices, batch_size=p
            )
            mean = []
            var = []
            for tensors in scdl:
                x = tensors[_CONSTANTS.X_KEY]
                y = tensors[_CONSTANTS.LABELS_KEY]
                out = self.module.inference(x, y)
                mean_, var_ = out["qz_m"], out["qz_v"]
                mean += [mean_.cpu()]
                var += [var_.cpu()]

            mean_vprior[ct], var_vprior[ct] = np.array(torch.cat(mean)), np.array(
                torch.cat(var)
            )

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
