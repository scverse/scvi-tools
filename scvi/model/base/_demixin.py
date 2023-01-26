import logging
from typing import Optional, Sequence, Union

import numpy as np
import torch
from anndata import AnnData
from sklearn.covariance import EllipticEnvelope
from torch.distributions import Categorical, Normal
from tqdm import tqdm

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
    module's (e.g., VAE) methods and associate signatures.
    """

    @torch.inference_mode()
    def get_normalized_expression_iw(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 25,
        n_samples_overall: int = None,
        batch_size: Optional[int] = 64,
        filter_cells: bool = True,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        return_numpy: Optional[bool] = False,
        n_mc_samples_px: int = 500,
        n_cells_per_chunk: Optional[int] = 32,
        truncation=False,
        # n_mc=3,
        max_chunks: Optional[int] = 10,
    ) -> np.ndarray:
        """
        Computes importance-weighted expression levels within a subpopulation.

        There are three majors steps to obtain the expression levels.
        A first optional step consists in filtering out outlier cells, using the cells' latent representation.
        Next, we infer how the data should be split to compute expression levels.
        In particular, if the considered subpopulation is big, we divide the indices into K subgroups of cells

        Finally, we compute importance-weighted expression levels within each chunk, for each subgroup,
        which we then concatenate and return.

        Parameters
        ----------
        adata
            Anndata to use, defaults to None, by default None
        indices
            Indices of the subpopulation, by default None
        n_samples
            Number of samples to use per cell.
        n_samples_overall
            Number of samples to use overall. Overrides n_samples if not None.
        batch_size
            Batch size of the data loader
        filter_cells
            Whether cells should be filtered using outlier detection, by default True
        transform_batch
            Batch to use, by default None
        return_numpy
             Whether numpy should be returned, by default False
        n_mc_samples_px
            Number of overall samples per cell used to compute the marginal likelihood
        n_cells_per_chunk
            Number of cells in each subgroup, by default 500
        max_chunks
            Maximum number of subgroups to use, by default None


        Returns
        -------
        res
            Importance weighted expression levels of shape (Number of samples, number of genes)
        """
        # Step 1: Determine effective indices to use
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
            indices_ = self.filter_outlier_cells(
                adata=adata,
                indices=indices_,
                batch_size=batch_size,
            )

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

        normalized_exprs = []
        for chunk in tqdm(cell_chunks):
            # for _ in range(n_mc):
            normalized_exprs.append(
                self._importance_weighted_expression(
                    adata=adata,
                    indices=chunk,
                    n_samples=n_samples_per_cell,
                    batch_size=batch_size,
                    n_mc_samples_px=n_mc_samples_px,
                    truncation=truncation,
                ).numpy()
            )
        normalized_exprs = np.concatenate(normalized_exprs, 0)
        idx = np.random.choice(
            np.arange(len(normalized_exprs)), size=n_samples_overall, replace=True
        )
        normalized_exprs = normalized_exprs[idx]
        return normalized_exprs

    @torch.no_grad()
    def _importance_weighted_expression(
        self,
        adata: AnnData,
        indices: Sequence,
        n_samples: int,
        n_mc_samples_px: int = 5000,
        batch_size: int = 64,
        truncation=False,
    ) -> dict:
        """
        Obtain gene expression and densities.

        Computes gene normalized expression samples, as well as
        variational posterior densities and likelihoods for each samples and each
        cell.

        For each cell of the dataset, we sample `n_samples` posterior samples.
        We then compute the likelihood of each sample and for each cell.
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
            h = generative_outputs["px"].scale
            n_genes = h.shape[-1]
            h = h.reshape(-1, n_genes).cpu()

            zs.append(inference_outputs["z"].reshape(-1, self.module.n_latent).cpu())
            qzs_m.append(inference_outputs["qz"].loc.cpu())
            qzs_std.append(inference_outputs["qz"].scale.cpu())
            hs.append(h)

            _log_px_zs = self.evaluate_likelihood_from_z(scdl, inference_outputs["z"])
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
            batch_size=batch_size,
            observation_specific=True,
            n_mc_samples_per_pass=250,
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
        sampled_particles = (
            Categorical(logits=log_probs.unsqueeze(0))
            .sample((ws.shape[0],))
            .squeeze(-1)
        )
        return hs[sampled_particles]

    @torch.no_grad()
    def evaluate_likelihood_from_z(self, scdl, z: torch.Tensor) -> torch.Tensor:
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
        z
            Latent codes for which to evaluate the likelihood.
        """
        z_reshaped = z.unsqueeze(1).reshape(-1, 1, self.module.n_latent)
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

    def filter_outlier_cells(
        self, adata: AnnData, indices: Sequence, batch_size: int
    ) -> Sequence:
        """
        Filter outlier cells indexed by indices.

        Parameters
        ----------
        adata:
            Anndata containing the observations.
        indices:
            Indices characterizing the considered subpopulation.
        batch_size:
            Batch-size to use to compute the latent representation.
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
