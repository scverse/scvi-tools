import logging
from typing import Optional

import numpy as np
import torch
from anndata import AnnData
from torch.utils.data import DataLoader, TensorDataset

from scvi import _CONSTANTS
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import TrainingPlan
from scvi.model.base import BaseModelClass, RNASeqMixin, VAEMixin
from scvi.module import VAEC

logger = logging.getLogger(__name__)


class CondSCVI(RNASeqMixin, VAEMixin, BaseModelClass):
    """
    Conditional version of single-cell Variational Inference, used for hierarchical deconvolution of spatial transcriptomics data.

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
    use_gpu
        Use the GPU or not.
    **model_kwargs
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
        n_latent: int = 10,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        use_gpu: bool = True,
        **module_kwargs,
    ):
        super(CondSCVI, self).__init__(adata, use_gpu=use_gpu)

        self.module = VAEC(
            n_input=self.summary_stats["n_vars"],
            n_labels=self.summary_stats["n_labels"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **module_kwargs,
        )
        self._model_summary_string = (
            "Conditional SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
        )
        self.init_params_ = self._get_init_params(locals())

    @property
    def _plan_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

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
            raise RuntimeError("Please train the model first.")

        adata = self._validate_anndata(adata)

        mean_vprior = np.zeros(
            (self.summary_stats["n_labels"], p, self.module.n_latent)
        )
        var_vprior = np.zeros((self.summary_stats["n_labels"], p, self.module.n_latent))
        key = adata.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["original_key"]
        mapping = adata.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"]
        for ct in range(self.summary_stats["n_labels"]):
            # pick p cells
            local_indices = np.random.choice(
                np.where(adata.obs[key] == mapping[ct])[0], p
            )
            # get mean and variance from posterior
            scdl = self._make_scvi_dl(adata=adata, indices=local_indices, batch_size=p)
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

    @torch.no_grad()
    def generate_from_latent(
        self,
        z: Optional[np.ndarray] = None,
        labels: Optional[np.ndarray] = None,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Return the scaled parameter of the NB for every cell.

        Parameters
        ----------
        z
            Numpy array with latent space
        labels
            Numpy array with labels
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        gene_expression
        """
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")

        dl = DataLoader(
            TensorDataset(torch.tensor(z), torch.tensor(labels, dtype=torch.long)),
            batch_size=128,
        )  # create your dataloader

        rate = []
        for tensors in dl:
            px_rate = self.module.generative(
                tensors[0], torch.ones((tensors[0].shape[0], 1)), tensors[1]
            )["px_scale"]
            rate += [px_rate.cpu()]
        return np.array(torch.cat(rate))
