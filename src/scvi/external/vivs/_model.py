from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import torch
import xarray as xr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField, ObsmField
from scvi.external.vivs._constants import VIVS_REGISTRY_KEYS
from scvi.external.vivs._module import VIVSModule
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin, VAEMixin
from scvi.module import VAE
from scvi.utils import dependencies, setup_anndata_dsp, track

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class VIVS(VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """VIVS: calibrated identification of feature dependencies in multiomics :cite:p:`BoyeauVIVS24`.

    Identifies which genes in ``X`` are conditionally dependent on an external response
    ``Y`` (e.g. protein expression, niche composition) using a conditional randomization
    test (CRT), with a deep generative model of ``X`` as the knockoff sampler.

    Parameters
    ----------
    adata
        AnnData object registered via :meth:`~scvi.external.VIVS.setup_anndata`.
    n_hidden
        Number of hidden units in the generative VAE (ignored if ``x_model`` is given).
    n_latent
        Latent dimensionality of the generative VAE (ignored if ``x_model`` is given).
    x_likelihood
        Gene-expression likelihood for the generative VAE: ``"nb"``, ``"zinb"``, or ``"poisson"``
        (ignored if ``x_model`` is given).
    xy_linear
        If ``True``, the importance-score net is a linear model instead of an MLP.
    xy_include_batch_in_input
        Whether to concatenate one-hot batch to the importance-score net's input.
    x_model
        An already-trained :class:`~scvi.model.SCVI` or :class:`~scvi.external.SCVIVA`
        instance (or any model whose ``.module`` is a :class:`~scvi.module.VAE` or
        subclass), registered on data compatible with ``adata``. When given, its trained
        module is reused (frozen) as the knockoff sampler, and VIVS's own generative
        training phase is skipped entirely. Models with a fundamentally different module
        API (e.g. :class:`~scvi.model.DestVI`, :class:`~scvi.external.RESOLVI`,
        :class:`~scvi.external.GIMVI`) are not supported here.
    **module_kwargs
        Additional keyword arguments passed to ``VIVSModule``.
    """

    _module_cls = VIVSModule

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        x_likelihood: Literal["nb", "zinb", "poisson"] = "nb",
        xy_linear: bool = False,
        xy_include_batch_in_input: bool = False,
        x_model: BaseModelClass | None = None,
        **module_kwargs,
    ):
        super().__init__(adata)
        summary_stats = self.summary_stats

        x_module = None
        if x_model is not None:
            if not x_model.is_trained:
                raise ValueError(
                    "`x_model` must already be trained (call `.train()` on it) before "
                    "being passed to VIVS."
                )
            if not isinstance(x_model.module, VAE):
                raise ValueError(
                    "`x_model.module` must be a `scvi.module.VAE` (or subclass, e.g. "
                    "SCVIVA's `nicheVAE`) so it can be reused as VIVS's knockoff sampler. "
                    f"Got `{type(x_model.module).__name__}` (from `{type(x_model).__name__}`), "
                    "whose inference()/generative() API is not compatible."
                )
            x_model_var_names = np.asarray(x_model.get_var_names())
            adata_var_names = np.asarray(self.get_var_names())
            if not np.array_equal(x_model_var_names, adata_var_names):
                raise ValueError(
                    "`x_model` was registered on data with different genes (or a different "
                    "gene order) than `adata`. VIVS reuses `x_model`'s encoder/decoder "
                    "weights column-for-column, so gene order must match exactly."
                )
            x_model_batch_mapping = x_model.adata_manager.get_state_registry(
                REGISTRY_KEYS.BATCH_KEY
            ).categorical_mapping
            adata_batch_mapping = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.BATCH_KEY
            ).categorical_mapping
            if not np.array_equal(x_model_batch_mapping, adata_batch_mapping):
                raise ValueError(
                    "`x_model` was registered with a different batch category mapping than "
                    "`adata`. Reusing its module would feed batch indices through the wrong "
                    "one-hot columns for the frozen encoder/decoder."
                )
            x_module = x_model.module

        self.module = self._module_cls(
            n_input=summary_stats.n_vars,
            n_responses=summary_stats.n_Y,
            n_batch=summary_stats.n_batch,
            x_model_kwargs={
                "n_hidden": n_hidden,
                "n_latent": n_latent,
                "gene_likelihood": x_likelihood,
            },
            xy_linear=xy_linear,
            xy_include_batch_in_input=xy_include_batch_in_input,
            x_module=x_module,
            **module_kwargs,
        )
        self._model_summary_string = "VIVS model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        y_obsm_key: str,
        layer: str | None = None,
        batch_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        y_obsm_key
            Key in ``adata.obsm`` for the response(s) ``Y`` whose conditional dependence on
            gene expression ``X`` is being tested (e.g. protein expression, niche composition).
        %(param_layer)s
        %(param_batch_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, None),
            ObsmField(VIVS_REGISTRY_KEYS.Y_KEY, y_obsm_key),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def train(
        self,
        max_epochs: int | None = None,
        x_max_epochs: int | None = None,
        xy_max_epochs: int | None = None,
        train_size: float = 0.9,
        validation_size: float | None = None,
        batch_size: int = 128,
        early_stopping: bool = False,
        **kwargs,
    ):
        """Train VIVS in two sequential phases.

        Phase 1 fits the generative VAE over ``X`` (skipped entirely if a pretrained
        ``x_model`` was supplied at construction). Phase 2 freezes it and fits the
        importance-score net for ``Y|X``. This order is required for CRT validity: the
        knockoff sampler must not be contaminated by information about ``Y``.
        """
        if not self.module.x_module_is_pretrained:
            self.module._phase = "x"
            super().train(
                max_epochs=x_max_epochs or max_epochs,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                early_stopping=early_stopping,
                **kwargs,
            )
            self.module.x_module.requires_grad_(False)
            self.module.x_module.eval()

        self.module._phase = "xy"
        super().train(
            max_epochs=xy_max_epochs or max_epochs,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            early_stopping=early_stopping,
            **kwargs,
        )
        self.is_trained_ = True

    @torch.inference_mode()
    def predict_t(
        self,
        adata: AnnData | None = None,
        indices=None,
        batch_size: int = 128,
    ) -> np.ndarray:
        """Raw per-cell importance-score-net predictions (no CRT knockoff perturbation)."""
        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        device = self.module.device
        results = []
        for tensors in dataloader:
            x = tensors[REGISTRY_KEYS.X_KEY].to(device)
            y = tensors[VIVS_REGISTRY_KEYS.Y_KEY].to(device)
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].to(device)
            xy_input = self.module.xy_input(x, batch_index)
            out = self.module.xy_module(xy_input, y)
            results.append(out["all_loss"].cpu().numpy())
        return np.concatenate(results, axis=0)

    @torch.inference_mode()
    def _encode_for_knockoffs(
        self, x: torch.Tensor, batch_index: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Run the frozen generative VAE's encoder once per batch.

        The original JAX implementation samples ``z`` a single time per batch and only
        resamples the decoder's noise across MC knockoff draws (`vivs/_vivs.py:227-228,302`
        on VIVS `main`: ``z_rng`` is split once, ``randomize(...)`` is called once per batch,
        outside the MC loop). Re-drawing ``z`` fresh on every MC sample (as an earlier draft
        of this helper did) folds extra encoder-uncertainty variance into every null draw and
        is NOT what the reference algorithm does — call this once per batch, then call
        `_sample_knockoffs` many times with the same `z`/`library`.
        """
        inference_out = self.module.inference(x=x, batch_index=batch_index)
        return inference_out["z"], inference_out["library"]

    @torch.inference_mode()
    def _sample_knockoffs(
        self, z: torch.Tensor, library: torch.Tensor, batch_index: torch.Tensor
    ) -> torch.Tensor:
        """Sample one conditional replacement of X from the frozen generative VAE's decoder.

        `z`/`library` must come from `_encode_for_knockoffs`, called ONCE per batch outside
        the MC loop — only the decoder's `px.sample()` varies across MC knockoff draws,
        matching the reference algorithm exactly (see `_encode_for_knockoffs`'s docstring).
        """
        generative_out = self.module.generative(z=z, library=library, batch_index=batch_index)
        return generative_out["px"].sample()

    @staticmethod
    def _crt_pvalue(
        obs_t: np.ndarray, tilde_t: np.ndarray, n_mc_samples: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """One-sided CRT p-value with BH correction, shared by every hypothesis-testing method.

        Parameters
        ----------
        obs_t
            Observed statistic, shape ``(..., n_responses)``.
        tilde_t
            Null statistics, shape ``(n_mc_samples, ..., n_responses)``.
        """
        pval = (1.0 + (obs_t >= tilde_t).sum(0)) / (1.0 + n_mc_samples)
        padj = np.stack(
            [multipletests(pval[..., r], method="fdr_bh")[1] for r in range(pval.shape[-1])],
            axis=-1,
        )
        return pval, padj

    @torch.inference_mode()
    def get_importance(
        self,
        adata: AnnData | None = None,
        indices=None,
        batch_size: int = 128,
        n_mc_samples: int = 500,
        use_vmap: Literal["auto", True, False] = "auto",
    ) -> dict:
        """Conditional-randomization-test importance of each gene for each response.

        Parameters
        ----------
        use_vmap
            Whether to vectorize the per-gene resampling loop with :func:`torch.vmap`.
            ``"auto"`` enables it when the number of genes is below 2000 (mirrors the
            original's own recommended gene-filtering ceiling). Disable if you hit an
            out-of-memory error.
        """
        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        n_genes = self.summary_stats.n_vars
        n_responses = self.summary_stats.n_Y
        use_vmap = use_vmap if use_vmap != "auto" else n_genes < 2000
        device = self.module.device

        obs_t_total = torch.zeros(n_responses, device=device)
        tilde_t_total = torch.zeros(n_mc_samples, n_genes, n_responses, device=device)
        n_obs = 0

        for tensors in dataloader:
            x = tensors[REGISTRY_KEYS.X_KEY].to(device)
            y = tensors[VIVS_REGISTRY_KEYS.Y_KEY].to(device)
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].to(device)
            batch_n = x.shape[0]
            n_obs += batch_n

            obs_all_loss = self.module.xy_module(self.module.xy_input(x, batch_index), y)[
                "all_loss"
            ]
            obs_t_total += obs_all_loss.sum(0)

            z, library = self._encode_for_knockoffs(x, batch_index)  # once per batch
            px_sample = self._sample_knockoffs(z, library, batch_index)  # (batch_n, n_genes)

            if use_vmap:
                try:
                    tilde_t_batch = self._get_importance_vmap(x, px_sample, batch_index, y)
                except RuntimeError as e:
                    raise RuntimeError(
                        "Out of memory while vmapping over genes. Try setting use_vmap=False."
                    ) from e
            else:
                tilde_t_batch = torch.stack(
                    [
                        self._compute_gene_statistic(x, px_sample[:, g], g, batch_index, y)
                        for g in range(n_genes)
                    ],
                    dim=0,
                )  # (n_genes, n_responses)
            # Only one MC sample per batch pass in this minimal loop-based version; repeat
            # `n_mc_samples` times by resampling. This mirrors `n_mc_per_pass=1` in the original.
            tilde_t_total[0] += tilde_t_batch
            for k in range(1, n_mc_samples):
                px_sample_k = self._sample_knockoffs(
                    z, library, batch_index
                )  # same z/library, fresh px draw
                if use_vmap:
                    try:
                        tilde_t_k = self._get_importance_vmap(x, px_sample_k, batch_index, y)
                    except RuntimeError as e:
                        raise RuntimeError(
                            "Out of memory while vmapping over genes. Try setting use_vmap=False."
                        ) from e
                else:
                    tilde_t_k = torch.stack(
                        [
                            self._compute_gene_statistic(x, px_sample_k[:, g], g, batch_index, y)
                            for g in range(n_genes)
                        ],
                        dim=0,
                    )
                tilde_t_total[k] += tilde_t_k

        obs_t = (obs_t_total / n_obs).cpu().numpy()
        null_t = (tilde_t_total / n_obs).cpu().numpy()
        pval, padj = self._crt_pvalue(obs_t, null_t, n_mc_samples)
        return {"obs_ts": obs_t, "null_ts": null_t, "pvalues": pval, "padj": padj}

    @torch.inference_mode()
    def get_cell_scores(
        self,
        gene_ids: list[int],
        response_ids: list[int] | None = None,
        adata: AnnData | None = None,
        indices=None,
        batch_size: int | None = None,
        n_mc_samples: int | None = None,
    ) -> dict:
        """Per-cell (unsummed) importance scores for a specific set of genes."""
        adata = self._validate_anndata(adata)
        batch_size = batch_size or 128
        n_mc_samples = n_mc_samples or 500
        response_ids = (
            response_ids if response_ids is not None else list(range(self.summary_stats.n_Y))
        )
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        device = self.module.device

        tilde_t_mean_chunks, obs_t_chunks = [], []
        for tensors in dataloader:
            x = tensors[REGISTRY_KEYS.X_KEY].to(device)
            y = tensors[VIVS_REGISTRY_KEYS.Y_KEY].to(device)
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].to(device)

            obs_all_loss = self.module.xy_module(self.module.xy_input(x, batch_index), y)[
                "all_loss"
            ][:, response_ids]
            obs_t = obs_all_loss  # (batch_n, n_responses_selected)

            z, library = self._encode_for_knockoffs(x, batch_index)  # once per batch
            tilde_t_sum = torch.zeros(x.shape[0], len(response_ids), device=device)
            for _ in range(n_mc_samples):
                px_sample = self._sample_knockoffs(
                    z, library, batch_index
                )  # fresh px draw, same z
                x_perturbed = x.clone()
                x_perturbed[..., gene_ids] = px_sample[..., gene_ids]
                xy_input = self.module.xy_input(x_perturbed, batch_index)
                all_loss = self.module.xy_module(xy_input, y)["all_loss"][:, response_ids]
                tilde_t_sum += all_loss
            tilde_t_mean_chunks.append((tilde_t_sum / n_mc_samples).cpu().numpy())
            obs_t_chunks.append(obs_t.cpu().numpy())

        return {
            "tilde_t_mean": np.concatenate(tilde_t_mean_chunks, axis=0),
            "obs_t": np.concatenate(obs_t_chunks, axis=0),
        }

    def _get_importance_vmap(
        self,
        x: torch.Tensor,
        x_tilde: torch.Tensor,
        batch_index: torch.Tensor,
        y: torch.Tensor,
    ) -> torch.Tensor:
        """Vectorized version of `_compute_gene_statistic`, mapped over the gene axis.

        Mirrors `scvi.external.mrvi.MRVI`'s vmap usage pattern: `torch.vmap` over one
        axis, `randomness="different"` since downstream computation involves stochastic
        submodules (dropout in the importance-score net). Uses a non-mutating
        `torch.where` substitution rather than in-place indexed assignment
        (`x_perturbed[..., gene_id] = ...`) — verified empirically that the in-place
        form raises `RuntimeError: vmap: index_put_(...)` as soon as a dropout/BatchNorm
        submodule is called downstream inside the vmapped function; `torch.where` does not.
        """
        n_genes = x.shape[-1]

        def _statistic_for_gene(gene_id: torch.Tensor, x_tilde_gene: torch.Tensor) -> torch.Tensor:
            mask = torch.nn.functional.one_hot(gene_id, n_genes).bool()
            x_tilde_expanded = x_tilde_gene.unsqueeze(-1).expand(-1, n_genes)
            x_perturbed = torch.where(mask, x_tilde_expanded, x)
            xy_input = self.module.xy_input(x_perturbed, batch_index)
            return self.module.xy_module(xy_input, y)["all_loss"].sum(0)

        gene_ids = torch.arange(n_genes, device=x.device)
        return torch.vmap(_statistic_for_gene, in_dims=(0, 1), randomness="different")(
            gene_ids, x_tilde
        )

    def _compute_gene_statistic(
        self,
        x: torch.Tensor,
        x_tilde_gene: torch.Tensor,
        gene_id: int,
        batch_index: torch.Tensor,
        y: torch.Tensor,
    ) -> torch.Tensor:
        """Substitute one gene with its knockoff, recompute xy statistic (summed over cells)."""
        x_perturbed = x.clone()
        x_perturbed[..., gene_id] = x_tilde_gene
        xy_input = self.module.xy_input(x_perturbed, batch_index)
        return self.module.xy_module(xy_input, y)["all_loss"].sum(0)

    @torch.inference_mode()
    def get_gene_correlations(
        self,
        adata: AnnData | None = None,
        indices=None,
        batch_size: int = 128,
    ) -> np.ndarray:
        """Gene-by-gene correlation matrix of the decoder's normalized expression scale."""
        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        n_genes = self.summary_stats.n_vars
        device = self.module.device

        x_sum = torch.zeros(n_genes, device=device)
        xx_sum = torch.zeros(n_genes, n_genes, device=device)
        n_obs = 0
        for tensors in dataloader:
            x = tensors[REGISTRY_KEYS.X_KEY].to(device)
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].to(device)
            inference_out = self.module.inference(x=x, batch_index=batch_index)
            generative_out = self.module.generative(
                z=inference_out["z"], library=inference_out["library"], batch_index=batch_index
            )
            scale = generative_out["px"].scale
            x_sum += scale.sum(0)
            xx_sum += torch.bmm(scale.unsqueeze(-1), scale.unsqueeze(1)).sum(0)
            n_obs += x.shape[0]

        x_mean = (x_sum / n_obs).unsqueeze(0)
        cov = xx_sum / n_obs - x_mean.T @ x_mean
        inv_std = 1.0 / torch.sqrt(torch.diag(cov))
        d = torch.diag(inv_std)
        corr = d @ cov @ d
        return corr.cpu().numpy()

    @dependencies("fastcluster")
    def get_gene_groupings(
        self,
        adata: AnnData | None = None,
        method: str = "complete",
        return_z=False,
        n_clusters_list: list[int] | None = None,
    ):
        """Hierarchically cluster genes by their decoder-scale correlation.

        Parameters
        ----------
        method
            Linkage method for hierarchical clustering.
        return_z
            Whether to also return the linkage matrix and computed gene order.
        n_clusters_list
            Cluster-count resolutions to compute a partition for.
        """
        import fastcluster

        assert n_clusters_list is not None
        adata = self._validate_anndata(adata)
        corr = self.get_gene_correlations(adata=adata)
        pseudo_dist = 1 - corr
        pseudo_dist = (pseudo_dist + pseudo_dist.T) / 2
        pseudo_dist = np.clip(pseudo_dist, a_min=0.0, a_max=100.0)
        pseudo_dist = pseudo_dist - np.diag(np.diag(pseudo_dist))
        dist_vec = squareform(pseudo_dist, checks=False)

        Z = fastcluster.linkage(dist_vec, method=method)
        Z = hierarchy.optimal_leaf_ordering(Z, dist_vec)
        gene_order = hierarchy.leaves_list(Z)
        gene_order = adata.var_names[gene_order].values

        n_genes = self.summary_stats.n_vars
        if not isinstance(n_clusters_list, list):
            n_clusters_list = [n_clusters_list]
        gene_groupings = []
        for n_cluster in n_clusters_list:
            if n_cluster >= n_genes:
                continue
            cluster_assignments = hierarchy.fcluster(Z, n_cluster, criterion="maxclust") - 1
            gene_groupings.append(cluster_assignments)

        if return_z:
            return gene_groupings, Z, gene_order
        return gene_groupings

    @torch.inference_mode()
    def get_hier_importance(
        self,
        n_clusters_list: list[int],
        adata: AnnData | None = None,
        indices=None,
        batch_size: int = 128,
        gene_groupings: list[np.ndarray] | None = None,
        gene_order=None,
        clustering_method: str = "complete",
        use_vmap: Literal["auto", True, False] = "auto",
        n_mc_samples: int = 500,
        silent: bool = False,
    ) -> xr.Dataset:
        """Hierarchical CRT: gene importance at multiple gene-group resolutions.

        First clusters genes by decoder-scale correlation at several resolutions (unless
        pre-computed groupings are given), then re-runs the CRT with group-level knockoff
        substitution at each resolution. See ``docs/user_guide/models/vivs.md`` for the
        full statistical description.

        Parameters
        ----------
        silent
            If ``True``, disables the progress bar tracking MC-sample/batch iterations.
        """
        adata = self._validate_anndata(adata)
        n_genes = self.summary_stats.n_vars
        n_responses = self.summary_stats.n_Y
        use_vmap = use_vmap if use_vmap != "auto" else n_genes < 2000
        device = self.module.device

        if gene_groupings is None:
            gene_groupings, _, gene_order = self.get_gene_groupings(
                adata=adata,
                n_clusters_list=n_clusters_list,
                return_z=True,
                method=clustering_method,
            )
        elif gene_order is None:
            gene_order = adata.var_names.values
        gene_groupings = [*gene_groupings, np.arange(n_genes).astype(np.int32)]
        gene_groups_oh = [
            torch.nn.functional.one_hot(torch.as_tensor(g).long()).float().to(device)
            for g in gene_groupings
        ]
        group_sizes = [oh.shape[1] for oh in gene_groups_oh]

        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        obs_t_total = torch.zeros(n_responses, device=device)
        tilde_t_totals = [
            torch.zeros(n_mc_samples, sz, n_responses, device=device) for sz in group_sizes
        ]
        n_obs = 0

        n_batches = len(dataloader)
        pbar = track(
            range(n_batches * n_mc_samples),
            description="Computing hierarchical importance",
            disable=silent,
        )
        pbar_iter = iter(pbar)

        for tensors in dataloader:
            x = tensors[REGISTRY_KEYS.X_KEY].to(device)
            y = tensors[VIVS_REGISTRY_KEYS.Y_KEY].to(device)
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY].to(device)
            n_obs += x.shape[0]

            obs_t_total += self.module.xy_module(self.module.xy_input(x, batch_index), y)[
                "all_loss"
            ].sum(0)

            z, library = self._encode_for_knockoffs(x, batch_index)  # once per batch
            for k in range(n_mc_samples):
                x_tilde = self._sample_knockoffs(z, library, batch_index)  # fresh px draw, same z
                for res_idx, group_oh in enumerate(gene_groups_oh):
                    stat = self._compute_group_statistics(
                        x, x_tilde, group_oh, batch_index, y, use_vmap
                    )
                    tilde_t_totals[res_idx][k] += stat
                next(pbar_iter, None)

        obs_t = (obs_t_total / n_obs).cpu().numpy()
        tilde_t = [(t / n_obs).cpu().numpy() for t in tilde_t_totals]
        return self._construct_hier_results(
            obs_t, tilde_t, gene_groupings, group_sizes, gene_order, adata
        )

    def _compute_group_statistics(
        self,
        x: torch.Tensor,
        x_tilde: torch.Tensor,
        group_oh: torch.Tensor,
        batch_index: torch.Tensor,
        y: torch.Tensor,
        use_vmap: bool,
    ) -> torch.Tensor:
        """Substitute each gene-group with its knockoff (soft mask) and recompute the statistic."""

        def _statistic_for_group(group_mask: torch.Tensor) -> torch.Tensor:
            x_perturbed = x * (1.0 - group_mask) + x_tilde * group_mask
            xy_input = self.module.xy_input(x_perturbed, batch_index)
            return self.module.xy_module(xy_input, y)["all_loss"].sum(0)

        if use_vmap:
            try:
                return torch.vmap(_statistic_for_group, in_dims=1, randomness="different")(
                    group_oh
                )
            except RuntimeError as e:
                raise RuntimeError(
                    "Out of memory while vmapping over gene groups. Try setting use_vmap=False."
                ) from e
        return torch.stack(
            [_statistic_for_group(group_oh[:, i]) for i in range(group_oh.shape[1])]
        )

    def _construct_hier_results(
        self,
        obs_t: np.ndarray,
        tilde_t: list[np.ndarray],
        gene_groupings: list[np.ndarray],
        group_sizes: list[int],
        gene_order,
        adata: AnnData,
    ) -> xr.Dataset:
        """Assemble the multi-resolution p-value/padj cube, reusing `_crt_pvalue`."""
        n_mc_samples = tilde_t[0].shape[0]
        response_names = np.arange(obs_t.shape[-1])
        datasets = []
        for res_idx, resolution in enumerate(group_sizes):
            pvals, padjs = self._crt_pvalue(obs_t, tilde_t[res_idx], n_mc_samples)
            gene_clusters = gene_groupings[res_idx]
            pvals_gene = pvals[gene_clusters]
            padjs_gene = padjs[gene_clusters]
            coords = {
                "gene_name": adata.var_names.values,
                "feature": response_names,
                "resolution": resolution,
                "resolution_idx": res_idx,
            }
            datasets.append(
                xr.Dataset(
                    {
                        "pval": (["gene_name", "feature"], pvals_gene),
                        "padj": (["gene_name", "feature"], padjs_gene),
                        "cluster_assignment": (["gene_name"], gene_clusters),
                    },
                    coords=coords,
                )
            )
        return xr.concat(datasets, dim="resolution").reindex(gene_name=gene_order)
