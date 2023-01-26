import logging
from typing import Optional, Sequence, Union

import numpy as np
import torch
from anndata import AnnData
from sklearn.covariance import EllipticEnvelope
from torch.distributions import Categorical, Normal

from scvi import REGISTRY_KEYS
from scvi.model._utils import _get_batch_code_from_category

logger = logging.getLogger(__name__)
Number = Union[int, float]


class DEMixin:
    """
    DE module relying on Importance-sampling.

    Mixin for using
    importance-weighted DE content.
    This however requires some additional structure on the
    module's (e.g., VAE) methods and associate signatures
    """

    @torch.no_grad()
    def get_subpopulation_expression(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 25,
        n_samples_overall: int = None,
        batch_size: Optional[int] = 64,
        filter_cells: bool = True,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        return_numpy: Optional[bool] = False,
        marginal_n_samples_per_pass: int = 250,
        n_mc_samples_px: int = 500,
        n_cells_per_chunk: Optional[int] = 500,
        truncation=False,
        n_mc=1,
        max_chunks: Optional[int] = None,
    ) -> np.ndarray:
        """
        Computes importance-weighted expression levels within a subpopulation.

        There are three majors steps to obtain the expression levels.
        A first optional step consists in filtering out outlier cells, using the cells' latent representation.
        Next, we infer how the data should be split to compute expression levels.
        In particular, if the considered subpopulation is big, we divide the indices
        into K chunks

        Finally, we compute importance-weighted expression levels within each chunk, in `_inference_loop`,
        which we then concatenate and return.

        Parameters
        ----------
        adata
            Anndata to use, defaults to None, by default None
        indices
            Indices of the subpopulation, by default None
        n_samples
            Indices of the subpopulation, by default 25
        n_samples_overall
            Indices of the subpopulation, by default None
        batch_size
            Batch size of the data loader, by default 64
        filter_cells
            Whether cells should be filtered using outlier detection, by default True
        transform_batch
            Batch to use, by default None
        return_numpy
             Whether numpy should be returned, by default False
        marginal_n_samples_per_pass
            Number of samples per pass to compute the marginal likelihood, by default 500
        n_mc_samples_px
            Number of overall samples per cell used to compute the marginal likelihood, by default 5000
        n_cells_per_chunk
            Number of cells to use in each minibatch, by default 500
        max_chunks
            Maximum number of chunks to use, by default None


        Returns
        -------
        res
            Importance weighted expression levels of shape (Number of samples, number of genes)
        """
        # Step 1: Determine effective indices to use
        # adata = self._validate_anndata(adata)
        if transform_batch is not None:
            adata_manager = self.get_anndata_manager(adata, required=True)
            transform_batch_val = _get_batch_code_from_category(
                adata_manager, transform_batch
            )[0]
            observed_batches = adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY)[
                indices
            ].squeeze()
            if indices.shape != observed_batches.shape:
                raise ValueError("Discrepancy between # of indices and # of batches")
            indices_ = indices[observed_batches == transform_batch_val]
        else:
            indices_ = indices
        if len(indices_) == 0:
            n_genes = adata.n_vars
            return np.array([]).reshape((0, n_genes))
        if filter_cells:
            indices_ = self._filter_cells(
                adata=adata,
                indices=indices_,
                batch_size=batch_size,
            )

        # Step 2.
        # Determine number of cell chunks
        # Because of the quadratic complexity we split cells in smaller chunks when
        # the cell population is too big
        # This ensures that the method remains scalable when looking at very large cell populations
        n_cells = indices_.shape[0]
        logger.debug(f"n cells {n_cells}")
        n_cell_chunks = int(np.ceil(n_cells / n_cells_per_chunk))
        np.random.seed(0)
        np.random.shuffle(indices_)
        cell_chunks = np.array_split(indices_, n_cell_chunks)[:max_chunks]
        n_cells_used = np.concatenate(cell_chunks).shape[0]
        # Determine number of samples to generate per cell
        if n_samples_overall is not None:
            n_samples_per_cell = int(1 + np.ceil(n_samples_overall / n_cells_used))
        else:
            n_samples_per_cell = n_samples
            n_samples_overall = n_samples_per_cell * n_cells_used
        n_samples_per_cell = np.minimum(n_samples_per_cell, 30)
        logger.debug(f"n samples per cell {n_samples_per_cell}")

        # Step 3
        res = []
        for chunk in cell_chunks:
            # for chunk in cell_indices:
            logger.debug(f"n cells in chunk {chunk.shape[0]}")
            for _ in range(n_mc):
                res.append(
                    self._importance_weighted_expression(
                        adata=adata,
                        indices=chunk,
                        n_samples=n_samples_per_cell,
                        batch_size=batch_size,
                        marginal_n_samples_per_pass=marginal_n_samples_per_pass,
                        n_mc_samples_px=n_mc_samples_px,
                        truncation=truncation,
                    ).numpy()
                )
            # logger.debug(res[-1].shape)
        res = np.concatenate(res, 0)
        idx = np.arange(len(res))
        idx = np.random.choice(idx, size=n_samples_overall, replace=True)
        res = res[idx]
        return res

    @torch.no_grad()
    def _importance_weighted_expression(
        self,
        adata: AnnData,
        indices: Sequence,
        n_samples: int,
        marginal_n_samples_per_pass: int = 500,
        n_mc_samples_px: int = 5000,
        batch_size: int = 64,
        truncation=False,
    ) -> dict:
        """
        Obtain gene expression and densities.

        Computes gene normalized expression samples, as well as
        variational posterior densities and likelihoods for each samples and each
        cell.

        The functioning of the method is a follows.
        For each cell of the dataset, we sample `n_samples` posterior samples.
        We then compute the likelihood of each sample and FOR EACH CELL, using
        `_evaluate_likelihood`.
        After concatenating likelihoods and samples, we derive overall
        importance weights and importance weighted expression levels.

        Parameters
        ----------
        adata
            Considered anndataset
        indices
            Indices of the subpopulation
        n_samples
            Number of posterior samples per cell
        marginal_n_samples_per_pass
            Number of samples per pass to compute the marginal likelihood, by default 500
        n_mc_samples_px
            Number of overall samples per cell used to compute the marginal likelihood, by default 5000
        batch_size
            Number of cells per minibatch, by default 64

        Returns
        -------
        dict
            Containing expression levels, z samples as well as associated densities
        """
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, shuffle=False
        )

        zs = []
        qzs_m = []
        qzs_std = []
        hs = []
        log_px_zs = []
        for tensors in scdl:
            inference_outputs, generative_outputs, = self.module.forward(
                tensors,
                inference_kwargs=dict(n_samples=n_samples),
                compute_loss=False,
            )
            z = inference_outputs["z"].reshape(-1, self.module.n_latent).cpu()
            h = generative_outputs["px"].scale
            n_genes = h.shape[-1]
            h = h.reshape(-1, n_genes).cpu()

            zs.append(z)
            qzs_m.append(inference_outputs["qz"].loc.cpu())
            qzs_std.append(inference_outputs["qz"].scale.cpu())
            hs.append(h)

            _log_px_zs = self._evaluate_likelihood(scdl, inference_outputs)
            log_px_zs.append(_log_px_zs)
        log_px_zs = torch.cat(log_px_zs, 0)

        zs = torch.cat(zs, dim=0)  # shape n_samples, n_cells, n_latent
        hs = torch.cat(hs, dim=0)
        qzs_m = torch.cat(qzs_m, dim=0)
        qzs_std = torch.cat(qzs_std, dim=0)

        _zs = zs.unsqueeze(1)  # shape (overall samples, 1, n_latent)
        log_qz = Normal(qzs_m, qzs_std).log_prob(_zs).sum(-1)
        log_pz = (
            Normal(torch.zeros_like(zs), torch.ones_like(zs))
            .log_prob(zs)
            .sum(-1)
            .squeeze()
        )

        log_px = self.get_marginal_ll(
            adata=adata,
            indices=indices,
            n_mc_samples=n_mc_samples_px,
            n_samples_per_pass=marginal_n_samples_per_pass,
            batch_size=batch_size,
            observation_specific=True,
        )

        importance_weight = torch.logsumexp(
            log_pz.unsqueeze(1)
            + log_px_zs
            - log_px
            - torch.logsumexp(log_qz, 1, keepdims=True),
            dim=1,
        )

        if truncation:
            tau = torch.logsumexp(importance_weight, 0) - np.log(
                importance_weight.shape[0]
            )
            importance_weight[importance_weight <= tau] = tau

        log_probs = importance_weight - torch.logsumexp(importance_weight, 0)
        ws = log_probs.exp()
        n_samples_overall = ws.shape[0]
        logger.debug(f"Max ratio {ws.max().item()}")
        windices = (
            Categorical(logits=log_probs.unsqueeze(0))
            .sample((n_samples_overall,))
            .squeeze(-1)
        )
        return hs[windices]

    @torch.no_grad()
    def _evaluate_likelihood(self, scdl, inference_outputs) -> torch.Tensor:
        r"""
        Derive required likelihoods.

        Computes :math:`p(x \mid z)`, :math:`q(z \mid x)` as well as :math:`p(z)` for
        each cell :math:`x` contained in `scdl` and predetermined
        posterior samples :math:`z` in `inference_outputs`.
        These quantities are necessary to evalute subpopulation-wide importance weights.

        Parameters
        ----------
        scdl
             Dataset containing cells of interest
        inference_outputs
            Output of module containing the latent variables z of interest

        Returns
        -------
        indices
            Likelihoods of shape (number of cells, number of posterior samples)
        """
        z_reshaped = (
            inference_outputs["z"].unsqueeze(1).reshape(-1, 1, self.module.n_latent)
        )
        n_samples_loop = z_reshaped.shape[0]
        log_px_zs = []
        for _tensors in scdl:
            inference_inputs = self.module._get_inference_input(_tensors)
            n_cells = inference_inputs["batch_index"].shape[0]
            zs = z_reshaped.expand(n_samples_loop, n_cells, self.module.n_latent)
            log_px_zs.append(
                self.module.estimate_likelihood(
                    tensors=_tensors, z=zs, library=None
                ).cpu()
            )
        log_px_zs = torch.cat(log_px_zs, 1)
        return log_px_zs

    def _filter_cells(
        self, adata: AnnData, indices: Sequence, batch_size: int
    ) -> Sequence:
        """
        Filter outlier cells indexed by indices.

        Parameters
        ----------
        adata
            Anndata containing the observations.
        indices
            Indices characterizing the considered subpopulation.
        batch_size
            Batch-size to use to compute the latent representation.

        Returns
        -------
        indices
            Indices to keep for differential expression
        """
        qz_m = self.get_latent_representation(
            adata=adata,
            indices=indices,
            give_mean=True,
            batch_size=batch_size,
        )
        if (qz_m.ndim != 2) or (qz_m.shape[0] != len(indices)):
            raise ValueError("Dimension mismatch of variational density means")
        try:
            idx_filt = EllipticEnvelope().fit_predict(qz_m)
            idx_filt = idx_filt == 1
        except ValueError:
            logger.warning("Could not properly estimate Cov!, using all samples")
            idx_filt = np.ones(qz_m.shape[0], dtype=bool)
        if (idx_filt == 1).sum() <= 1:
            idx_filt = np.ones(qz_m.shape[0], dtype=bool)
        try:
            logger.debug(
                "n cells total {}, n cells filtered {}".format(
                    idx_filt.shape[0], idx_filt.sum()
                )
            )
            indices = indices[idx_filt]
        except IndexError:
            raise IndexError(idx_filt)
        return indices
