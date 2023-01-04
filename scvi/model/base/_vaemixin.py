import logging
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import torch
from anndata import AnnData

from scvi.utils import unsupported_in_latent_mode

from ._log_likelihood import compute_elbo, compute_reconstruction_error

logger = logging.getLogger(__name__)


class VAEMixin:
    """Univseral VAE methods."""

    @torch.inference_mode()
    @unsupported_in_latent_mode
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
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        elbo = compute_elbo(self.module, scdl)
        return -elbo

    @torch.inference_mode()
    @unsupported_in_latent_mode
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
        if indices is None:
            indices = np.arange(adata.n_obs)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        if hasattr(self.module, "marginal_ll"):
            log_lkl = 0
            for tensors in scdl:
                log_lkl += self.module.marginal_ll(tensors, n_mc_samples=n_mc_samples)
        else:
            raise NotImplementedError(
                "marginal_ll is not implemented for current model. "
                "Please raise an issue on github if you need it."
            )
        n_samples = len(indices)
        return log_lkl / n_samples

    @torch.inference_mode()
    @unsupported_in_latent_mode
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
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        reconstruction_error = compute_reconstruction_error(self.module, scdl)
        return reconstruction_error

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: Optional[int] = None,
        return_dist: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
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
        return_dist
            Return the distribution parameters of the latent variables rather than their sampled values.
            If `True`, ignores `give_mean` and `mc_samples`.

        Returns
        -------
        Low-dimensional representation for each cell or a tuple containing its mean and variance.
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        latent = []
        latent_qzm = []
        latent_qzv = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            if "qz" in outputs:
                qz = outputs["qz"]
            else:
                qz_m, qz_v = outputs["qz_m"], outputs["qz_v"]
                qz = torch.distributions.Normal(qz_m, qz_v.sqrt())
            if give_mean:
                # does each model need to have this latent distribution param?
                if self.module.latent_distribution == "ln":
                    samples = qz.sample([mc_samples])
                    z = torch.nn.functional.softmax(samples, dim=-1)
                    z = z.mean(dim=0)
                else:
                    z = qz.loc
            else:
                z = outputs["z"]

            latent += [z.cpu()]
            latent_qzm += [qz.loc.cpu()]
            latent_qzv += [qz.scale.square().cpu()]
        return (
            (torch.cat(latent_qzm).numpy(), torch.cat(latent_qzv).numpy())
            if return_dist
            else torch.cat(latent).numpy()
        )
