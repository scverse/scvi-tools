from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import torch

from scvi.data._utils import _validate_adata_dataloader_input
from scvi.utils import unsupported_if_adata_minified

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    import numpy.typing as npt
    from anndata import AnnData
    from torch import Tensor
    from torch.distributions import Distribution


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
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        return_mean: bool = True,
        **kwargs,
    ) -> float:
        """Compute the evidence lower bound (ELBO) on the data.

        The ELBO is the reconstruction error plus the Kullback-Leibler (KL) divergences between the
        variational distributions and the priors. It is different from the marginal log-likelihood;
        specifically, it is a lower bound on the marginal log-likelihood plus a term that is
        constant with respect to the variational distribution. It still gives good insights on the
        modeling of the data and is fast to compute.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
            Ignored if ``dataloader`` is not ``None``.
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        return_mean
            Whether to return the mean of the ELBO or the ELBO for each observation.
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

        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )
        else:
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )

        return -compute_elbo(self.module, dataloader, return_mean=return_mean, **kwargs)

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
            Ignored if ``dataloader`` is not ``None``.
        n_mc_samples
            Number of Monte Carlo samples to use for the estimator. Passed into the module's
            ``marginal_ll`` method.
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``.
        return_mean
            Whether to return the mean of the marginal log-likelihood or the marginal-log
            likelihood for each observation.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        **kwargs
            Additional keyword arguments to pass into the module's ``marginal_ll`` method.

        Returns
        -------
        If ``True``, returns the mean marginal log-likelihood. Otherwise returns a tensor of shape
        ``(n_obs,)`` with the marginal log-likelihood for each observation.

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
        else:
            _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )
        else:
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
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
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        return_mean: bool = True,
        **kwargs,
    ) -> dict[str, float]:
        r"""Compute the reconstruction error on the data.

        The reconstruction error is the negative log likelihood of the data given the latent
        variables. It is different from the marginal log-likelihood, but still gives good insights
        on the modeling of the data and is fast to compute. This is typically written as
        :math:`p(x \mid z)`, the likelihood term given one posterior sample.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
            Ignored if ``dataloader`` is not ``None``
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        return_mean
            Whether to return the mean reconstruction loss or the reconstruction loss
            for each observation.
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

        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )
        else:
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )

        return compute_reconstruction_error(
            self.module, dataloader, return_mean=return_mean, **kwargs
        )

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
            Ignored if ``dataloader`` is not ``None``
        give_mean
            If ``True``, returns the mean of the latent distribution. If ``False``, returns an
            estimate of the mean using ``mc_samples`` Monte Carlo samples.
        mc_samples
            Number of Monte Carlo samples to use for the estimator for distributions with no
            closed-form mean (e.g., the logistic normal distribution). Not used if ``give_mean`` is
            ``True`` or if ``return_dist`` is ``True``.
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``
        return_dist
            If ``True``, returns the mean and variance of the latent distribution. Otherwise,
            returns the mean of the latent distribution.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.

        Returns
        -------
        An array of shape ``(n_obs, n_latent)`` if ``return_dist`` is ``False``. Otherwise, returns
        a tuple of arrays ``(n_obs, n_latent)`` with the mean and variance of the latent
        distribution.
        """
        from torch.distributions import Normal
        from torch.nn.functional import softmax

        from scvi.module._constants import MODULE_KEYS

        self._check_if_trained(warn=False)
        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )
        else:
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
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
