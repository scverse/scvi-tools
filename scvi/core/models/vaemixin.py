import logging
from typing import Optional, Sequence

import numpy as np
import torch
from anndata import AnnData

from scvi import _CONSTANTS
from scvi.core.trainers import UnsupervisedTrainer

logger = logging.getLogger(__name__)


class VAEMixin:
    def train(
        self,
        n_epochs: Optional[int] = None,
        train_size: float = 0.9,
        test_size: Optional[float] = None,
        lr: float = 1e-3,
        n_epochs_kl_warmup: int = 400,
        n_iter_kl_warmup: Optional[int] = None,
        frequency: Optional[int] = None,
        train_fun_kwargs: dict = {},
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        n_epochs
            Number of passes through the dataset.
        train_size
            Size of training set in the range [0.0, 1.0].
        test_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + test_size < 1`, the remaining cells belong to a validation set.
        lr
            Learning rate for optimization.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        n_iter_kl_warmup
            Number of minibatches for scaling term on KL divergence to go from 0 to 1.
            To use, set to not `None` and set `n_epochs_kl_warmup` to `None`.
        frequency
            Frequency with which metrics are computed on the data for train/test/val sets.
        train_fun_kwargs
            Keyword args for the train method of :class:`~scvi.core.trainers.UnsupervisedTrainer`.
        **kwargs
            Other keyword args for :class:`~scvi.core.trainers.UnsupervisedTrainer`.
        """
        train_fun_kwargs = dict(train_fun_kwargs)
        if self.is_trained_ is False:
            self.trainer = UnsupervisedTrainer(
                self.model,
                self.adata,
                train_size=train_size,
                test_size=test_size,
                n_iter_kl_warmup=n_iter_kl_warmup,
                n_epochs_kl_warmup=n_epochs_kl_warmup,
                frequency=frequency,
                use_cuda=self.use_cuda,
                **kwargs,
            )
            self.train_indices_ = self.trainer.train_set.indices
            self.test_indices_ = self.trainer.test_set.indices
            self.validation_indices_ = self.trainer.validation_set.indices
            self.history_ = self.trainer.history
        # for autotune
        if "n_epochs" not in train_fun_kwargs:
            if n_epochs is None:
                n_cells = self.adata.n_obs
                n_epochs = np.min([round((20000 / n_cells) * 400), 400])
            train_fun_kwargs["n_epochs"] = n_epochs
        if "lr" not in train_fun_kwargs:
            train_fun_kwargs["lr"] = lr

        logger.info("Training for {} epochs".format(n_epochs))
        self.trainer.train(**train_fun_kwargs)
        self.is_trained_ = True

    @torch.no_grad()
    def get_elbo(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> float:
        """
        Return the ELBO for the data.

        The ELBO is a lower bound on the log likelihood of the data used for optimization
        of VAEs. Note, this is not the negative ELBO, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        return -scdl.elbo()

    @torch.no_grad()
    def get_marginal_ll(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_mc_samples: int = 1000,
        batch_size: Optional[int] = None,
    ) -> float:
        """
        Return the marginal LL for the data.

        The computation here is a biased estimator of the marginal log likelihood of the data.
        Note, this is not the negative log likelihood, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_mc_samples
            Number of Monte Carlo samples to use for marginal LL estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        return -scdl.marginal_ll(n_mc_samples=n_mc_samples)

    @torch.no_grad()
    def get_reconstruction_error(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> float:
        r"""
        Return the reconstruction error for the data.

        This is typically written as :math:`p(x \mid z)`, the likelihood term given one posterior sample.
        Note, this is not the negative likelihood, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        return -scdl.reconstruction_error()

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Return the latent representation for each cell.

        This is denoted as :math:`z_n` in our manuscripts.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell
        """
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")

        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        for tensors in scdl:
            x = tensors[_CONSTANTS.X_KEY]
            b = tensors[_CONSTANTS.BATCH_KEY]
            z = self.model.sample_from_posterior_z(
                x, b, give_mean=give_mean, n_samples=mc_samples
            )
            latent += [z.cpu()]
        return np.array(torch.cat(latent))
