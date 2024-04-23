from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence

import numpy.typing as npt
import torch
from anndata import AnnData
from torch import Tensor

from scvi.utils import unsupported_if_adata_minified

logger = logging.getLogger(__name__)


class VAEMixin:
    """Universal variational auto-encoder (VAE) methods."""

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_elbo(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
    ) -> float:
        """Compute the evidence lower bound (ELBO) on the data.

        The ELBO is a lower bound on the log-likelihood of the data used for optimization of VAEs.

        The ELBO is a lower bound on the log likelihood of the data used for optimization
        of VAEs. Note, this is not the negative ELBO, higher is better.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
        batch_size
            Minibatch size for the forward pass. Only used if ``dataloader`` is ``None``. If
            ``None``, defaults to ``scvi.settings.batch_size``.
        dataloader
            An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
            minibatch formatted as expected by the model.
        **kwargs
            Additional keyword arguments to pass into the forward method of the module.

        Returns
        -------
        Evidence lower bound (ELBO) of the data.

        Notes
        -----
        This is not the negative ELBO, so higher is better.
        """
        from scvi.model.base._log_likelihood import compute_elbo

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return -compute_elbo(self.module, dataloader)

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_marginal_ll(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        n_mc_samples: int = 1_000,
        batch_size: int | None = None,
        return_mean: bool = True,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
        **kwargs,
    ) -> float | Tensor:
        """Compute the marginal log-likehood of the data.

        The computation here is a biased estimator of the marginal log-likelihood of the data.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
        n_mc_samples
            Number of Monte Carlo samples to use for the estimator.
        batch_size
            Minibatch size for the forward pass. Only used if ``dataloader`` is ``None``. If
            ``None``, defaults to ``scvi.settings.batch_size``.
        return_mean
            Whether to return the mean of the marginal log-likelihood or the marginal-log
            likelihood for each observation.
        dataloader
            An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
            minibatch formatted as expected by the model.
        **kwargs
            Additional keyword arguments to pass into the module's ``marginal_ll`` method.

        Returns
        -------
        A tensor of shape ``(n_obs,)`` with the marginal log-likelihood for each observation if
        ``return_mean`` is ``False``. Otherwise, returns the mean marginal log-likelihood.

        Notes
        -----
        This is not the negative log-likelihood, so higher is better.
        """
        from numpy import mean

        if not hasattr(self.module, "marginal_ll"):
            raise NotImplementedError(
                "The model's module must implement `marginal_ll` to compute the marginal "
                "log-likelihood."
            )

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        log_likelihoods: list[float | Tensor] = [
            self.module.marginal_ll(
                tensors, n_mc_samples=n_mc_samples, return_mean=return_mean, **kwargs
            )
            for tensors in dataloader
        ]

        if return_mean:
            return mean(log_likelihoods)
        else:
            return torch.cat(log_likelihoods, dim=0)

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_reconstruction_error(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
        **kwargs,
    ) -> float:
        r"""Compute the reconstruction error on the data.

        This is typically written as :math:`p(x \mid z)`, the likelihood term given one posterior
        sample. Note, this is not the negative likelihood, higher is better.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
        batch_size
            Minibatch size for the forward pass. Only used if ``dataloader`` is ``None``. If
            ``None``, defaults to ``scvi.settings.batch_size``.
        dataloader
            An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
            minibatch formatted as expected by the model.
        **kwargs
            Additional keyword arguments to pass into the forward method of the module.

        Returns
        -------
        Reconstruction error for the data.

        Notes
        -----
        This is not the negative reconstruction error, so higher is better.
        """
        from scvi.model.base._log_likelihood import compute_reconstruction_error

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return compute_reconstruction_error(self.module, dataloader, **kwargs)

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        mc_samples: int = 5_000,
        batch_size: int | None = None,
        return_dist: bool = False,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
    ) -> npt.NDArray | tuple[npt.NDArray, npt.NDArray]:
        """Compute the latent representation of the data.

        This is typically denoted as :math:`z_n`.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
        give_mean
            If ``True``, returns the mean of the latent distribution. If ``False``, returns an
            estimate of the mean using ``mc_samples`` Monte Carlo samples.
        mc_samples
            Number of Monte Carlo samples to use for the estimator for distributions with no
            closed-form mean (e.g., the logistic normal distribution). Not used if ``give_mean`` is
            ``True`` or if ``return_dist`` is ``True``.
        batch_size
            Minibatch size for the forward pass. Only used if ``dataloader`` is ``None``. If
            ``None``, defaults to ``scvi.settings.batch_size``.
        return_dist
            If ``True``, returns the mean and variance of the latent distribution. Otherwise,
            returns the mean of the latent distribution.
        dataloader
            An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
            minibatch formatted as expected by the model.

        Returns
        -------
        An array of shape ``(n_obs, n_latent)`` if ``return_dist`` is ``False``. Otherwise, returns
        a tuple of arrays ``(n_obs, n_latent)`` with the mean and variance of the latent
        distribution.
        """
        from torch.distributions import Distribution, Normal
        from torch.nn.functional import softmax

        from scvi.module._constants import MODULE_KEYS

        self._check_if_trained(warn=False)
        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        zs: list[Tensor] = []
        qz_means: list[Tensor] = []
        qz_vars: list[Tensor] = []
        for tensors in dataloader:
            outputs: dict[str, Tensor | Distribution | None] = self.module.inference(
                **self.module._get_inference_input(tensors)
            )

            if MODULE_KEYS.QZ_KEY in outputs:
                qz: Distribution = outputs.get(MODULE_KEYS.QZ_KEY)
                qzm: Tensor = qz.loc
                qzv: Tensor = qz.scale.square()
            else:
                qzm: Tensor = outputs.get(MODULE_KEYS.QZM_KEY)
                qzv: Tensor = outputs.get(MODULE_KEYS.QZV_KEY)
                qz: Distribution = Normal(qzm, qzv.sqrt())

            if return_dist:
                qz_means.append(qzm.cpu())
                qz_vars.append(qzv.cpu())
                continue

            z: Tensor = qzm if give_mean else outputs.get(MODULE_KEYS.Z_KEY)

            if give_mean and getattr(self.module, "latent_distribution", None) == "ln":
                samples = qz.sample([mc_samples])
                z = softmax(samples, dim=-1).mean(dim=0)

            zs.append(z.cpu())

        if return_dist:
            return torch.cat(qz_means).numpy(), torch.cat(qz_vars).numpy()
        else:
            return torch.cat(zs).numpy()
