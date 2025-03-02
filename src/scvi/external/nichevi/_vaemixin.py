from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import torch

from scvi.model.base import (
    VAEMixin,
)
from scvi.utils import unsupported_if_adata_minified

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from anndata import AnnData
    from torch import Tensor

logger = logging.getLogger(__name__)


class NicheVAEMixin(VAEMixin):
    """Adding niche losses computation to universal variational auto-encoder (VAE) methods."""

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_composition_error(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
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
        from ._log_likelihood import compute_composition_error

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return compute_composition_error(
            self.module, dataloader, return_mean=return_mean, **kwargs
        )

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_niche_error(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
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
        from ._log_likelihood import compute_niche_error

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return compute_niche_error(self.module, dataloader, return_mean=return_mean, **kwargs)
