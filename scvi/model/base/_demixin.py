from functools import partial
from typing import Dict, Iterable, Optional, Sequence, Union
import logging

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from torch.distributions import Categorical, Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.model._utils import (
    _get_batch_code_from_category,
    _get_var_names_from_setup_anndata,
    scrna_raw_counts_properties,
)
from scvi.model.base._utils import _de_core
from scvi.utils import track


logger = logging.getLogger(__name__)


class DEMixin:
    """Universal implementation for using
    importance-weighted DE content.
    This however requires some additional structure on the
    module's (e.g., VAE) methods and associate signatures
    """
    def lvm_de(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = None,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        eps: float = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        **kwargs
            Keyword args for :func:`scvi.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        col_names = _get_var_names_from_setup_anndata(adata)
        model_fn = partial(
            self.get_population_expression,
            return_numpy=True,
            batch_size=batch_size,
        )
        result = _de_core(
            adata,
            model_fn,
            groupby,
            group1,
            group2,
            idx1,
            idx2,
            all_stats,
            scrna_raw_counts_properties,
            col_names,
            mode,
            batchid1,
            batchid2,
            delta,
            batch_correction,
            fdr_target,
            silent,
            eps=eps,
            **kwargs,
        )

        return result



    @torch.no_grad()
    def get_population_expression(
        self,
        adata=None,
        indices=None,
        n_samples: int = 25,
        n_samples_overall: int = None,
        batch_size: int = 64,
        transform_batch: Optional[str] = None,
        return_numpy: Optional[bool] = False,
    ):
        """Returns scales and latent variable in a given subpopulation characterized by `indices`
        using importance sampling

        :param adata: Anndata to use, defaults to None
        :param indices: Indices of the subpopulation, defaults to None
        :param n_samples: Number of posterior samples per cell, defaults to 25
        :param n_samples_overall: Number of posterior samples to use in total, defaults to None.
        If provided, this parameter overrides `n_samples`.
        :param batch_size: Batch size of the data loader, defaults to 64
        :param transform_batch: Batch to use, defaults to None
        :return_numpy: Whether numpy should be returned
        """
        adata = self._validate_anndata(adata)
        if transform_batch is not None:
            adata_key = self.scvi_setup_dict_["data_registry"]["batch_indices"][
                "attr_key"
            ]
            observed_batches = adata[indices].obs[adata_key].values
            transform_batch_val = _get_batch_code_from_category(adata, transform_batch)[
                0
            ]
            indices_ = indices[observed_batches == transform_batch_val]
        else:
            indices_ = indices
        scdl = self._make_data_loader(
            adata=adata, indices=indices_, batch_size=batch_size
        )
        n_cells = scdl.indices.shape[0]

        if n_samples_overall is not None:
            n_samples_per_cell = int(np.ceil(n_samples_overall / n_cells))
        else:
            n_samples_per_cell = n_samples
            n_samples_overall = n_samples_per_cell * n_cells

        inference_outs = self._inference_loop(scdl, n_samples_per_cell)
        log_px_zs = inference_outs["log_px_zs"]
        log_qz = inference_outs["log_qz"]
        log_pz = inference_outs["log_pz"]
        hs = inference_outs["hs"]

        log_px = self.get_marginal_ll(
            adata=adata,
            indices=indices_,
            n_mc_samples=5000,
        )
        log_mixture = torch.logsumexp(log_px_zs - log_px, dim=-1)
        log_Q = torch.logsumexp(log_qz, dim=-1)
        log_target = log_pz + log_mixture - log_Q
        importance_weight = log_target - log_Q
        log_probs = importance_weight - torch.logsumexp(importance_weight, 0)
        logger.debug("ESS: {}".format(1 / (log_probs **2).sum().item()))
        indices = (
            Categorical(logits=log_probs.unsqueeze(0))
            .sample((n_samples_overall,))
            .squeeze(-1)
        )
        res = hs[indices]
        if return_numpy:
            res = res.numpy()
        return res

    @torch.no_grad()
    def _inference_loop(
        self,
        scdl,
        n_samples: int,
        # transform_batch: str = None
    ):
        """Returns population wide unweighted posterior samples, as well as
        variational posterior densities and likelihoods for each samples and each
        cell.

        :param scdl: Dataloader over subpopulation of interest
        :param n_samples: Number of samples per cell
        :param transform_batch: Batch of use, defaults to None
        """
        n_cells = scdl.indices.shape[0]
        zs = []
        qzs_m = []
        qzs_v = []
        hs = []
        log_px_zs = []
        for tensors in track(scdl, description="Iterating over cells ..."):
            inference_outputs, generative_outputs, = self.module.forward(
                tensors,
                inference_kwargs=dict(n_samples=n_samples, return_densities=True),
                compute_loss=False,
            )
            z = inference_outputs["z"].cpu()
            zs.append(z)
            qzs_m.append(inference_outputs["qz_m"].cpu())
            qzs_v.append(inference_outputs["qz_v"].cpu())
            hs.append(generative_outputs["px_scale"].cpu())

            _log_px_zs = self._evaluate_likelihood(scdl, inference_outputs)
            log_px_zs.append(_log_px_zs)
        log_px_zs = torch.cat(log_px_zs, 0)
        zs = torch.cat(zs, dim=1)  # shape n_samples, n_cells, n_latent
        zs = zs.reshape(n_cells * n_samples, self.module.n_latent)

        hs = torch.cat(hs, dim=1)
        n_genes = hs.shape[-1]
        hs = hs.reshape(n_cells * n_samples, n_genes)
        qzs_m = torch.cat(qzs_m, dim=0)
        qzs_v = torch.cat(qzs_v, dim=0)

        _zs = zs.unsqueeze(1)  # shape (overall samples, 1, n_latent)
        log_qz = Normal(qzs_m, qzs_v.sqrt()).log_prob(_zs).sum(-1)
        log_pz = (
            Normal(torch.zeros_like(zs), torch.ones_like(zs))
            .log_prob(zs)
            .sum(-1)
            .squeeze()
        )
        return dict(
            log_px_zs=log_px_zs,
            log_qz=log_qz,
            log_pz=log_pz,
            zs=zs,
            hs=hs,
            qzs_m=qzs_m,
            qzs_v=qzs_v,
        )

    @torch.no_grad()
    def _evaluate_likelihood(self, scdl, inference_outputs):
        """Computes p(x \mid z), q(z \mid x) as well as p(z) for
        each cell x contained in `scdl` and predetermined
        posterior samples $z$ in `inference_outputs`.
        These quantities are necessary to evalute subpopulation-wide importance weights

        :param scdl: Dataset containing cells of interest
        :param inference_outputs: Inference outputs containing the considered samples
        """
        z_samples = inference_outputs["z"]
        if self.module.use_observed_lib_size:
            lib_key = "library"
        else:
            lib_key = "ql_m"
        _z = z_samples.unsqueeze(1).reshape(-1, 1, self.module.n_latent)
        _n_samples_loop = _z.shape[0]
        _log_px_zs = []
        for _tensors in scdl:
            # This is simply used to get a good library value for the cells we are looking at
            _inf_inputs = self.module._get_inference_input(_tensors)
            _n_cells = _inf_inputs["batch_index"].shape[0]

            __z = _z.expand(_n_samples_loop, _n_cells, self.module.n_latent)

            point_library = (
                self.module.inference(**_inf_inputs, return_densities=True)[lib_key]
                .squeeze(0)
                .expand(_n_samples_loop, _n_cells, 1)
            )

            inference_outputs["z"] = __z.to()
            inference_outputs["library"] = point_library.to()
            _log_px_zs.append(
                self.module.generative_evaluate(
                    tensors=_tensors, inference_outputs=inference_outputs
                )["log_px_latents"].cpu()
            )
        _log_px_zs = torch.cat(_log_px_zs, 1)
        return _log_px_zs
