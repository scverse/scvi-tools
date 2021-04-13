"""May be easier to do a double loop to avoid working with complex tensors like in the reproducibility codebase
"""

import logging
from typing import Optional

import torch
import torch.nn
import numpy as np
from torch.distributions import Normal, Categorical
from anndata import AnnData

from scvi._compat import Literal
from scvi.external.wscvi._module import WVAE
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.utils import track
from scvi.model._utils import _get_batch_code_from_category


logger = logging.getLogger(__name__)


class WSCVI(
    RNASeqMixin, VAEMixin, ArchesMixin, UnsupervisedTrainingMixin, BaseModelClass
):
    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.0,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super(WSCVI, self).__init__(adata)
        self.n_genes = adata.shape[-1]
        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else None
        )
        self.module = WVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self._model_summary_string = (
            "SCVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            gene_likelihood,
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

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
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        n_cells = scdl.indices.shape[0]

        if n_samples_overall is not None:
            n_samples_per_cell = int(np.ceil(n_samples_overall / n_cells))
        else:
            n_samples_per_cell = n_samples
            n_samples_overall = n_samples_per_cell * n_cells

        inference_outs = self._inference_loop(scdl, n_samples_per_cell, transform_batch)
        log_px_zs = inference_outs["log_px_zs"]
        log_qz = inference_outs["log_qz"]
        log_pz = inference_outs["log_pz"]
        zs = inference_outs["zs"]
        hs = inference_outs["hs"]

        log_px = self._get_marginal(scdl)
        log_mixture = torch.logsumexp(log_px_zs - log_px, dim=-1)
        log_Q = torch.logsumexp(log_qz, dim=-1)
        log_target = log_pz + log_mixture - log_Q
        importance_weight = log_target - log_Q
        log_probs = importance_weight - torch.logsumexp(importance_weight, 0)
        indices = (
            Categorical(logits=log_probs.unsqueeze(0))
            .sample((n_samples_overall,))
            .squeeze(-1)
        )
        # return dict(
        #     h=hs[indices],
        #     z=zs[indices],
        # )
        res = hs[indices]
        if return_numpy:
            res = res.numpy()
        return res

    def _inference_loop(self, scdl, n_samples: int, transform_batch: str = None):
        """Returns population wide unweighted posterior samples, as well as
        variational posterior densities and likelihoods for each samples and each
        cell.

        :param scdl: Dataloader over subpopulation of interest
        :param n_samples: Number of samples per cell
        :param transform_batch: Batch of use, defaults to None
        """
        transform_batch = _get_batch_code_from_category(self.adata, transform_batch)[0]

        n_cells = scdl.indices.shape[0]
        if self.module.use_observed_lib_size:
            lib_key = "library"
        else:
            lib_key = "ql_m"
        zs = []
        qzs_m = []
        qzs_v = []
        hs = []
        log_px_zs = []
        for tensors in track(scdl, description="Iterating over cells ..."):
            inference_outputs, generative_outputs, = self.module.forward(
                tensors,
                inference_kwargs=dict(n_samples=n_samples),
                generative_kwargs=dict(transform_batch=transform_batch),
                compute_loss=False,
            )
            z = inference_outputs["z"].cpu()
            zs.append(z)
            qzs_m.append(inference_outputs["qz_m"].cpu())
            qzs_v.append(inference_outputs["qz_v"].cpu())
            hs.append(generative_outputs["px_scale"].cpu())

            assert inference_outputs[lib_key].ndim == 2
            _log_px_zs = self._evaluate_likelihood(scdl, inference_outputs)
            log_px_zs.append(_log_px_zs)
        log_px_zs = torch.cat(log_px_zs, 0)
        zs = torch.cat(zs, dim=1)  # shape n_samples, n_cells, n_latent
        zs = zs.reshape(n_cells * n_samples, self.module.n_latent)
        hs = torch.cat(hs, dim=1).reshape(n_cells * n_samples, self.n_genes)
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

    def _evaluate_likelihood(self, scdl, inference_outputs):
        """Computes p(x \mid z), q(z \mid x) as well as p(z) for
        each cell x contained in `scdl` and predetermined
        posterior samples $z$ in `inference_outputs`.
        These quantities are necessary to evalute subpopulation-wide importance weights

        :param scdl: Dataset containing cells of interest
        :param inference_outputs: Inference outputs containing the considered samples
        """
        z_samples = inference_outputs["z"]
        _z = z_samples.unsqueeze(1).reshape(-1, 1, self.module.n_latent)
        _n_samples_loop = _z.shape[0]
        _log_px_zs = []
        for _tensors in scdl:
            # This is simply used to get a good library value for the cells we are looking at
            _inf_inputs = self.module._get_inference_input(_tensors)
            _n_cells = _inf_inputs["batch_index"].shape[0]

            __z = _z.expand(_n_samples_loop, _n_cells, self.module.n_latent)
            point_library = (
                self.module.inference(**_inf_inputs)["point_library"]
                .squeeze(0)
                .expand(_n_samples_loop, _n_cells, 1)
            )

            inference_outputs["z"] = __z
            inference_outputs["library"] = point_library
            # _z =
            _gen_tensors = self.module._get_generative_input(
                _tensors, inference_outputs
            )
            _log_px_zs.append(
                self.module.generative(**_gen_tensors)["log_px_latents"].cpu()
            )
        _log_px_zs = torch.cat(_log_px_zs, 1)
        return _log_px_zs

    @torch.no_grad()
    def _get_marginal(
        self, scdl, n_mc_samples=5000, n_samples_per_pass=100, transform_batch=None
    ):
        """Computes cell-specific log evidence for each cell in `scdl`

        :param scdl: Dataloader to use
        :param n_mc_samples: Number of Monte Carlo samples, defaults to 5000
        :param n_samples_per_pass: Number of samples per forward pass, defaults to 100
        :param transform_batch: Batch to use, defaults to None
        """
        log_px = []
        n_iterations = int(n_mc_samples // n_samples_per_pass)
        for tensors in scdl:
            log_pxi = []
            for _ in range(n_iterations):
                inference_outputs, generative_outputs, = self.module.forward(
                    tensors,
                    inference_kwargs=dict(n_samples=n_samples_per_pass),
                    generative_kwargs=dict(transform_batch=transform_batch),
                    compute_loss=False,
                )
                log_joint = (
                    generative_outputs["log_px_latents"]
                    + generative_outputs["log_pz"]
                    + generative_outputs["log_pl"]
                )
                log_q = inference_outputs["log_ql"] + inference_outputs["log_qz"]
                log_ratios = log_joint - log_q
                log_pxi.append(log_ratios.cpu())
            log_pxi = torch.cat(log_pxi, 0)
            log_px.append(log_pxi)
        log_px_ = torch.cat(log_px, dim=1)
        log_px_ = torch.logsumexp(log_px_, 0)
        return log_px_
