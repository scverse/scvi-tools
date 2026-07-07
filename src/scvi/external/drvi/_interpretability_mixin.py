from __future__ import annotations

import inspect
import itertools
import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from torch.nn import functional as F

import scvi
from scvi import REGISTRY_KEYS
from scvi.external.drvi._constants import DRVI_MODULE_KEYS
from scvi.module._constants import MODULE_KEYS

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any, Literal

    from anndata import AnnData
    from torch import Tensor

logger = logging.getLogger(__name__)


def _flatten_to_1d(arr: np.ndarray | sparse.csr_matrix) -> np.ndarray:
    """Flatten a dense or sparse array to a 1D dense array."""
    dense = arr.todense() if sparse.issparse(arr) else arr
    return np.asarray(dense).flatten()


def _sparse_std(x: sparse.csr_matrix, axis: int = 0, ddof: int = 0) -> np.ndarray:
    """Standard deviation of a sparse matrix without densifying."""
    mean_sq = np.asarray(x.power(2).mean(axis=axis)).squeeze(axis=axis)
    sq_mean = np.asarray(np.power(x.mean(axis=axis), 2)).squeeze(axis=axis)
    var = mean_sq - sq_mean
    if ddof > 0:
        n = x.shape[axis]
        var = var * (n / (n - ddof))
    return np.sqrt(np.maximum(var, 0))


class InterpretabilityMixin:
    """Interpretability analyses for DRVI's latent factors (splits).

    Quantifies how much each latent split contributes to the
    reconstructed expression, both in-distribution (over the data) and out-of-distribution (by
    traversing each latent dimension), and outputs the per-gene scores for downstream use.
    Relies on the additive split decoder exposing per-split parameters in ``inspect_mode`` (see
    :class:`~scvi.external.drvi.DecoderDRVI`) and on the :class:`GenerativeMixin` iterators.
    """

    def _get_n_cats_per_cov(self) -> list[int]:
        """Number of categories for each categorical covariate (empty if none registered)."""
        if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry:
            state = self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY)
            return list(state.n_cats_per_key)
        return []

    def _get_norm_of_splits(self, latent: Tensor) -> Tensor:
        """Reduce latent ``(..., n_latent)`` to L2 norms ``(..., n_split)``."""
        n_split = self.module.n_split_latent
        n_latent = latent.shape[-1]
        return latent.reshape(*latent.shape[:-1], n_split, n_latent // n_split).norm(dim=-1)

    @torch.inference_mode()
    def iterate_on_effect_of_splits_within_distribution(
        self,
        adata: AnnData | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        add_to_counts: float = 1.0,
        deterministic: bool = True,
        directional: bool = True,
        **kwargs: Any,
    ):
        """Yield ``(effect_tensor, latent)`` per minibatch: each split's effect on every gene.

        Low-level building block behind the in-distribution scores: it runs the autoencoder over
        the data and converts the per-split decoder parameters into a per-gene effect. The effect
        formula depends on ``module.split_aggregation`` (``"logsumexp"`` or ``"mean"``).

        Parameters
        ----------
        adata
            AnnData to run on. Defaults to the model's registered AnnData.
        dataloader
            Custom minibatch iterator (e.g. from an out-of-core datamodule) used instead of
            building one from ``adata``. Exactly one of ``adata`` / ``dataloader`` may be given.
        add_to_counts
            Pseudo-count (relative to a ``1e6``-count cell) used by the ``"logsumexp"`` effect for
            numerical stability.
        deterministic
            If ``True`` (default), use the posterior mean for the bottleneck (no sampling).
        directional
            If ``True`` (default), split each effect into the factor's positive and negative
            directions. Requires ``n_split_latent == n_latent``.
        kwargs
            Forwarded to :meth:`iterate_on_ae_output` (e.g. ``indices``, ``batch_size``).

        Yields
        ------
        effect_tensor
            Per-split, per-gene effect, ``(n_cells, n_splits, n_genes)`` (or
            ``(n_cells, 2, n_splits, n_genes)`` when ``directional``).
        latent
            Posterior-mean latent for the minibatch, ``(n_cells, n_latent)``.
        """
        if directional and self.module.n_latent != self.module.n_split_latent:
            raise NotImplementedError(
                "Directional in-distribution interpretability requires one split per latent."
            )

        for inference_outputs, generative_outputs in self.iterate_on_ae_output(
            adata=adata,
            dataloader=dataloader,
            deterministic=deterministic,
            **kwargs,
        ):
            # posterior mean, n_cells x n_latent
            latent = inference_outputs[MODULE_KEYS.QZ_KEY].loc

            if self.module.split_aggregation == "logsumexp":
                # softmax(x) == softmax(x + c); the library / -log(K) constants cancel, but the
                # add_to_counts offset depends on the total decoder effect (a 1e6-count cell).
                log_mean_params = generative_outputs[
                    DRVI_MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY
                ]  # n_samples x n_splits x n_genes
                total_effect_per_cell = torch.logsumexp(log_mean_params, dim=[-2, -1])
                log_add_to_counts = total_effect_per_cell + np.log(add_to_counts / 1e6)
                log_mean_params = torch.cat(
                    [
                        log_mean_params,
                        log_add_to_counts.reshape(-1, 1, 1).expand(
                            -1, -1, log_mean_params.shape[-1]
                        ),
                    ],
                    dim=-2,
                )  # n_samples x (n_splits + 1) x n_genes
                effect_tensor = -torch.log(1 - F.softmax(log_mean_params, dim=-2)[:, :-1, :])
            elif self.module.split_aggregation == "mean":
                effect_tensor = torch.abs(
                    generative_outputs[DRVI_MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY]
                )
            else:
                raise NotImplementedError(
                    "Only 'logsumexp' and 'mean' aggregations are supported for interpretability."
                )

            if directional:
                pos_mask = (latent > 0).float().unsqueeze(-1)
                neg_mask = (latent < 0).float().unsqueeze(-1)
                effect_tensor = torch.stack(
                    [effect_tensor * pos_mask, effect_tensor * neg_mask],
                    dim=1,
                )  # n_samples x 2 x n_splits x n_genes

            yield effect_tensor, latent

    @torch.inference_mode()
    def get_reconstruction_effect_of_each_split(
        self,
        adata: AnnData | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        add_to_counts: float = 1.0,
        aggregate_over_cells: bool = True,
        deterministic: bool = True,
        directional: bool = False,
        **kwargs: Any,
    ) -> np.ndarray:
        """Effect of each split on reconstruction (summed over genes).

        Used by :meth:`set_latent_dimension_stats` to rank latent factors by how much they drive
        the reconstruction.

        Parameters
        ----------
        adata
            AnnData to run on. Defaults to the model's registered AnnData.
        dataloader
            Custom minibatch iterator used instead of building one from ``adata``. Exactly one of
            ``adata`` / ``dataloader`` may be given.
        add_to_counts
            Pseudo-count forwarded to the effect computation.
        aggregate_over_cells
            If ``True`` (default), sum the effect over all cells; otherwise return per-cell values.
        deterministic
            If ``True`` (default), use the posterior mean for the bottleneck (no sampling).
        directional
            If ``True``, split into the factor's positive and negative directions.
        kwargs
            Forwarded to :meth:`iterate_on_effect_of_splits_within_distribution`.

        Returns
        -------
        ``(n_splits,)`` (or ``(2, n_splits)`` if directional) when aggregating over cells, else
        ``(n_cells, n_splits)`` (or ``(n_cells, 2, n_splits)``).
        """
        store = None if aggregate_over_cells else []
        for effect_tensor, _latent in self.iterate_on_effect_of_splits_within_distribution(
            adata=adata,
            dataloader=dataloader,
            add_to_counts=add_to_counts,
            deterministic=deterministic,
            directional=directional,
            **kwargs,
        ):
            effect_share = effect_tensor.sum(dim=-1).detach().cpu()
            if aggregate_over_cells:
                effect_share = effect_share.sum(dim=0).to_dense()
                store = effect_share if store is None else store + effect_share
            else:
                store.append(effect_share)

        if aggregate_over_cells:
            return store.numpy(force=True)
        return torch.cat(store, dim=0).numpy(force=True)

    @torch.inference_mode()
    def set_latent_dimension_stats(
        self,
        embed: AnnData,
        adata: AnnData | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        vanished_threshold: float = 0.5,
    ) -> None:
        """Annotate ``embed.var`` with per-dimension reconstruction effect, ordering and stats.

        Typically the first interpretability step: it ranks the latent factors and records the
        per-dimension activation stats. When a split groups several latent dimensions
        (``n_split_latent < n_latent``), each split's reconstruction effect is broadcast to the
        latent dimensions it contains.

        Parameters
        ----------
        embed
            AnnData of the latent space (one ``var`` per latent factor, ``X`` holding the latent
            values). Annotated in place.
        adata
            AnnData to compute the reconstruction effect on. Defaults to the model's registered
            AnnData.
        dataloader
            Custom minibatch iterator used instead of building one from ``adata``. Exactly one of
            ``adata`` / ``dataloader`` may be given.
        vanished_threshold
            A dimension (or direction) is flagged as "vanished" (inactive) when its maximum
            absolute latent value stays below this threshold.

        Notes
        -----
        Adds the columns ``reconstruction_effect``, ``order``, ``max_value``, ``mean``, ``min``,
        ``max``, ``std``, ``std_abs``, ``title`` and ``vanished`` (+ per-direction vanished flags)
        to ``embed.var``.
        """
        n_latent = self.module.n_latent
        n_split = self.module.n_split_latent
        if "original_dim_id" not in embed.var:
            embed.var["original_dim_id"] = np.arange(embed.var.shape[0])

        # reconstruction effect is computed per split; when a split groups several latent dims
        # (n_split_latent < n_latent) broadcast each split's effect to the dims it contains.
        recon_per_split = self.get_reconstruction_effect_of_each_split(
            adata=adata, dataloader=dataloader
        )
        recon_per_dim = np.repeat(recon_per_split, n_latent // n_split)
        embed.var["reconstruction_effect"] = 0.0
        embed.var.loc[embed.var.sort_values("original_dim_id").index, "reconstruction_effect"] = (
            recon_per_dim
        )
        embed.var["order"] = (-embed.var["reconstruction_effect"]).argsort().argsort()

        embed.var["max_value"] = _flatten_to_1d(np.abs(embed.X).max(axis=0))
        embed.var["mean"] = _flatten_to_1d(embed.X.mean(axis=0))
        embed.var["min"] = _flatten_to_1d(embed.X.min(axis=0))
        embed.var["max"] = _flatten_to_1d(embed.X.max(axis=0))
        if sparse.issparse(embed.X):
            embed.var["std"] = _sparse_std(embed.X, axis=0)
            embed.var["std_abs"] = _sparse_std(np.abs(embed.X), axis=0)
        else:
            embed.var["std"] = embed.X.std(axis=0)
            embed.var["std_abs"] = np.abs(embed.X).std(axis=0)

        embed.var["title"] = "DR " + (1 + embed.var["order"]).astype(str)
        embed.var["vanished"] = embed.var["max_value"] < vanished_threshold
        embed.var["vanished_positive_direction"] = embed.var["max"] < vanished_threshold
        embed.var["vanished_negative_direction"] = embed.var["min"] > -vanished_threshold

    @torch.inference_mode()
    def get_effect_of_splits_within_distribution(
        self,
        adata: AnnData | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        add_to_counts: float = 1.0,
        deterministic: bool = True,
        directional: bool = True,
        aggregations: Sequence[Literal["max", "linear_weighted_mean", "exp_weighted_mean"]]
        | str = "ALL",
        skip_threshold: float = 1.0,
        **kwargs: Any,
    ) -> dict[str, np.ndarray]:
        """In-distribution per-split, per-gene effect scores.

        Aggregates each factor's per-gene effect over the observed data, weighting cells by how
        strongly the factor is activated.

        Parameters
        ----------
        adata
            AnnData to run on. Defaults to the model's registered AnnData.
        dataloader
            Custom minibatch iterator used instead of building one from ``adata``. Exactly one of
            ``adata`` / ``dataloader`` may be given.
        add_to_counts
            Pseudo-count forwarded to the effect computation.
        deterministic
            If ``True`` (default), use the posterior mean for the bottleneck (no sampling).
        directional
            If ``True`` (default), score the factor's positive and negative directions separately.
            Requires ``n_split_latent == n_latent``.
        aggregations
            Which cell-aggregation(s) to compute: ``"max"`` (peak effect across cells),
            ``"linear_weighted_mean"`` and ``"exp_weighted_mean"`` (activation-weighted means), or
            ``"ALL"`` (default) for all three. A single name or a sequence is accepted.
        skip_threshold
            Cells whose latent activation is below this value contribute (near) zero weight to the
            weighted-mean aggregations, so weakly-activated cells are effectively ignored.
        kwargs
            Forwarded to :meth:`iterate_on_effect_of_splits_within_distribution`.

        Returns
        -------
        A dict mapping each requested aggregation to a ``(n_splits, n_genes)`` array (or
        ``(2, n_splits, n_genes)`` when ``directional``).
        """
        if aggregations == "ALL":
            aggregations = ["max", "linear_weighted_mean", "exp_weighted_mean"]
        elif isinstance(aggregations, str):
            aggregations = [aggregations]

        store = {}
        for i, (effect_tensor, latent) in enumerate(
            self.iterate_on_effect_of_splits_within_distribution(
                adata=adata,
                dataloader=dataloader,
                add_to_counts=add_to_counts,
                deterministic=deterministic,
                directional=directional,
                **kwargs,
            )
        ):
            for aggregation in aggregations:
                if aggregation == "max":
                    effect_agg = (
                        effect_tensor.to_dense().amax(dim=0).detach().cpu().numpy(force=True)
                    )
                    if i == 0:
                        store[aggregation] = {"result": effect_agg}
                    else:
                        store[aggregation]["result"] = np.maximum(
                            store[aggregation]["result"], effect_agg
                        )
                elif aggregation in ("linear_weighted_mean", "exp_weighted_mean"):
                    if directional:
                        weights = torch.stack(
                            [
                                latent.clamp(min=skip_threshold) - skip_threshold,
                                (-latent).clamp(min=skip_threshold) - skip_threshold,
                            ],
                            dim=1,
                        ).unsqueeze(-1)  # n_samples x 2 x n_latent x 1
                    else:
                        # weighting based on norm of each split (abs when n_split == n_latent).
                        weights = self._get_norm_of_splits(latent).unsqueeze(-1)
                        weights = weights.clamp(min=skip_threshold) - skip_threshold
                    if aggregation == "exp_weighted_mean":
                        weights = torch.exp(weights) - 1.0
                    effect_agg = (
                        (effect_tensor * weights)
                        .sum(dim=0)
                        .to_dense()
                        .detach()
                        .cpu()
                        .numpy(force=True)
                    )
                    sum_weights = weights.sum(dim=0).detach().cpu().numpy(force=True)
                    if i == 0:
                        store[aggregation] = {"result": effect_agg, "sum_weights": sum_weights}
                    else:
                        store[aggregation]["result"] = store[aggregation]["result"] + effect_agg
                        store[aggregation]["sum_weights"] = (
                            store[aggregation]["sum_weights"] + sum_weights
                        )
                else:
                    raise NotImplementedError()

        results = {}
        for aggregation in aggregations:
            if aggregation == "max":
                results[aggregation] = store[aggregation]["result"]
            else:
                results[aggregation] = store[aggregation]["result"] / np.maximum(
                    store[aggregation]["sum_weights"], 1.0
                )
        return results

    @torch.inference_mode()
    def get_effect_of_splits_out_of_distribution(
        self,
        embed: AnnData,
        n_steps: int = 20,
        n_samples: int = 100,
        add_to_counts: float = 1.0,
        directional: bool = True,
        batch_size: int = scvi.settings.batch_size,
        seed: int | None = None,
    ) -> dict[str, np.ndarray]:
        """Out-of-distribution per-split, per-gene effect scores via latent traversal.

        Leverages DRVI's additive decoder: each latent factor ``Z_i`` is decoded by its own
        subnetwork ``f_i``, and the per-gene effects are aggregated over factors (``logsumexp``
        here). The effect of perturbing ``Z_i`` on a gene therefore depends only on ``f_i``, which
        is probed by traversing the dimension over discrete steps between its observed ``min``
        and ``max`` (from :meth:`set_latent_dimension_stats`), averaged over randomly sampled
        batch/covariate combinations. ``f_i``'s nonlinearity is summarized by its maximum
        contributions over the traversal.

        For each factor and gene, three log-fold-change (LFC) effect scores are returned, each
        ``(n_splits, n_genes)`` (or ``(2, n_splits, n_genes)`` when ``directional``, splitting the
        factor's negative ``[min, 0]`` and positive ``[0, max]`` traversal directions):

        * ``max_possible`` -- the **largest achievable** LFC: the factor's extreme contribution vs.
          its own minimum, with every other factor held at its minimum contribution.
          Overall shows how much a gene can be moved by the factor, regardless of other factors.
        * ``min_possible`` -- the same LFC but with every other factor held at its maximum
          contribution. Shrinks when a gene is driven by several factors.
          Overall specific is a factor's effect.
        * ``combined`` -- the product ``max_possible * min_possible``; large only when the factor
        is both strong and specific.

        ``add_to_counts`` is a pseudo-count (relative to a ``1e6``-count reference cell) folded
        into LFC calculations for numerical stability. See the DRVI paper :cite:p:`Moinfar2024`
        (Supplementary Note "Interpretation of Latent Factors via Additive Decoder Architecture")
        for the full derivation.

        This function Requires ``embed.var`` to contain ``min`` and ``max``.

        The batch/covariate combinations are sampled randomly. Pass ``seed`` for a reproducible
        sample; ``None`` (default) leaves the RNG untouched, so it still honors a previously set
        ``scvi.settings.seed``.
        """
        if self.module.n_latent != self.module.n_split_latent:
            raise NotImplementedError(
                "out-of-distribution interpretability requires one split per latent."
            )

        if n_steps % 2 != 0:
            raise ValueError(f"n_steps must be even, got {n_steps}.")
        dim_mins = np.minimum(embed.var["min"].values, 0.0)
        dim_maxs = np.maximum(embed.var["max"].values, 0.0)

        n_cat_total = [self.summary_stats.n_batch] + self._get_n_cats_per_cov()
        all_cat_combinations = np.asarray(
            list(itertools.product(*[range(n) for n in n_cat_total]))
        )
        # uses numpy's global RNG, so it is reproducible via ``scvi.settings.seed``; ``seed`` (when
        # given) reseeds it here, while ``seed=None`` leaves the RNG state untouched.
        if seed is not None:
            np.random.seed(seed)
        all_cat_combinations = np.random.permutation(all_cat_combinations)[:n_samples]

        n_combined = 0
        store = {"min_possible": None, "max_possible": None, "combined": None}
        for all_cats in all_cat_combinations:
            batch_val = all_cats[0].tolist()
            cat_vals = all_cats[1:].tolist()

            steps = np.linspace(1, 0, num=int(n_steps / 2), endpoint=False)
            span_values = np.concatenate(
                [steps[:, None] * dim_maxs[None, :], steps[:, None] * dim_mins[None, :]],
                axis=0,
            ).astype(np.float32)  # n_steps x n_latent

            effect_tensors = []
            for gen_output in self.iterate_on_decoded_latent_samples(
                z=span_values,
                lib=np.ones(n_steps) * 1e4,
                batch_values=np.full(n_steps, batch_val),
                cat_values=np.array([cat_vals] * n_steps) if cat_vals else None,
                cont_values=None,
                batch_size=batch_size,
                map_cat_values=False,
            ):
                effect_tensors.append(gen_output[DRVI_MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY])
            effect_tensors = torch.cat(effect_tensors, dim=0)  # n_steps x n_splits x n_genes
            effect_tensors = effect_tensors.reshape(
                2, int(n_steps / 2), effect_tensors.shape[1], effect_tensors.shape[2]
            )

            # extremes of each factor's log-space contribution f_i over the traversal.
            min_per_split = effect_tensors.amin(dim=[0, 1])  # n_splits x n_genes
            max_per_split = effect_tensors.amax(dim=[0, 1])  # n_splits x n_genes
            # max over the +/- half-ranges separately when directional, else over the full range
            directional_max_per_split = (
                effect_tensors.amax(dim=1) if directional else max_per_split
            )

            # add_to_counts pseudo-count (vs a 1e6-count cell) for numerical stability.
            log_add_to_counts = torch.logsumexp(max_per_split, dim=[0, 1]) + np.log(
                add_to_counts / 1e6
            )
            add_min = log_add_to_counts.reshape(1, 1).expand(1, min_per_split.shape[1])
            add_max = log_add_to_counts.reshape(1, 1).expand(1, max_per_split.shape[1])
            # All factors at maximum
            lse_min = torch.logsumexp(
                torch.cat([min_per_split, add_min], dim=0), dim=0, keepdim=True
            )
            # All factors at minimum
            lse_max = torch.logsumexp(
                torch.cat([max_per_split, add_max], dim=0), dim=0, keepdim=True
            )
            # effect_max-lfc: other factors held at their minimum (lse_min)
            max_possible = (
                torch.log(
                    torch.exp(lse_min)
                    - torch.exp(min_per_split)
                    + torch.exp(directional_max_per_split)
                )
                - lse_min
            )
            # effect_min-lfc: other factors held at their maximum (lse_max)
            min_possible = torch.log(
                torch.exp(lse_max)
                - torch.exp(max_per_split)
                + torch.exp(directional_max_per_split)
            ) - torch.log(torch.exp(lse_max) - torch.exp(max_per_split) + torch.exp(min_per_split))
            combined = max_possible * min_possible  # effect_combined

            if n_combined == 0:
                store["min_possible"] = min_possible
                store["max_possible"] = max_possible
                store["combined"] = combined
            else:
                store["min_possible"] = store["min_possible"] + min_possible
                store["max_possible"] = store["max_possible"] + max_possible
                store["combined"] = store["combined"] + combined
            n_combined += 1

        store = {k: (v / n_combined).detach().cpu().numpy() for k, v in store.items()}
        return store

    def calculate_interpretability_scores(
        self,
        embed: AnnData,
        methods: Sequence[str] | str = "OOD",
        directional: bool = True,
        add_to_counts: float = 1.0,
        inplace: bool = True,
        **kwargs: Any,
    ) -> dict[str, np.ndarray] | None:
        """Compute per-factor, per-gene interpretability scores and (optionally) store them.

        Top-level entry point for DRVI's latent-factor interpretability. It runs the in- and/or
        out-of-distribution analyses and collects their per-aggregation score matrices (one value
        per latent factor and gene) under a single, consistently named set of keys.

        Parameters
        ----------
        embed
            AnnData of the latent space (one ``var`` per latent factor), as annotated by
            :meth:`set_latent_dimension_stats`. Scores are written to its ``varm`` when
            ``inplace``.
        methods
            Which analyses to run: ``"IND"`` (in-distribution, over the data via
            :meth:`get_effect_of_splits_within_distribution`), ``"OOD"`` (default; out-of-
            distribution latent traversal via :meth:`get_effect_of_splits_out_of_distribution`), or
            ``"ALL"`` for both. A sequence of these is also accepted.
        directional
            If ``True`` (default), each factor's positive and negative directions are scored
            separately and the resulting keys gain a ``_positive`` / ``_negative`` suffix. Requires
            ``n_split_latent == n_latent``.
        add_to_counts
            Pseudo-count forwarded to the underlying effect computations (see those methods).
        inplace
            If ``True`` (default), store each score matrix in ``embed.varm`` and return ``None``;
            otherwise return them as a dict and leave ``embed`` untouched.
        kwargs
            Extra keyword arguments forwarded to the underlying method(s), each filtered to the
            arguments that method actually accepts. Useful ones include ``adata`` / ``dataloader``,
            ``deterministic`` and ``skip_threshold`` (IND), and ``n_steps``, ``n_samples``,
            ``seed`` and ``batch_size`` (OOD).

        Returns
        -------
        ``None`` when ``inplace=True`` (scores written to ``embed.varm``), otherwise a dict mapping
        each key to its ``(n_factors, n_genes)`` score matrix.

        Notes
        -----
        Keys are ``"{method}_{aggregation}[_positive|_negative]"``. The aggregations are ``max``,
        ``linear_weighted_mean`` and ``exp_weighted_mean`` for IND, and ``min_possible``,
        ``max_possible`` and ``combined`` for OOD (e.g. ``"OOD_combined_positive"``). These are the
        keys :meth:`get_interpretability_scores` reads back (default ``"OOD_combined"``).
        """
        if methods == "ALL":
            methods = ["IND", "OOD"]
        elif isinstance(methods, str):
            methods = [methods]
        calculate_ind = any(m.startswith("IND") for m in methods)
        calculate_ood = any(m.startswith("OOD") for m in methods)

        all_results = {}
        if calculate_ind:
            sig = inspect.signature(self.get_effect_of_splits_within_distribution)
            valid_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
            result_dict = self.get_effect_of_splits_within_distribution(
                directional=directional,
                add_to_counts=add_to_counts,
                aggregations="ALL",
                **valid_kwargs,
            )
            for key, value in result_dict.items():
                if directional:
                    all_results[f"IND_{key}_positive"] = value[0]
                    all_results[f"IND_{key}_negative"] = value[1]
                else:
                    all_results[f"IND_{key}"] = value
        if calculate_ood:
            sig = inspect.signature(self.get_effect_of_splits_out_of_distribution)
            valid_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
            result_dict = self.get_effect_of_splits_out_of_distribution(
                embed=embed, directional=directional, add_to_counts=add_to_counts, **valid_kwargs
            )
            for key, value in result_dict.items():
                if directional:
                    all_results[f"OOD_{key}_positive"] = value[0]
                    all_results[f"OOD_{key}_negative"] = value[1]
                else:
                    all_results[f"OOD_{key}"] = value

        if inplace:
            for key, value in all_results.items():
                if key in embed.varm:
                    logger.warning("Key %s already exists in embed.varm, overwriting.", key)
                embed.varm[key] = value
            return None
        return all_results

    def get_interpretability_scores(
        self,
        embed: AnnData,
        adata: AnnData,
        key: str = "OOD_combined",
        directional: bool = True,
        gene_symbols: str | None = None,
        order_col: str = "order",
        title_col: str = "title",
        hide_vanished: bool = True,
    ) -> pd.DataFrame:
        """Return interpretability scores as a genes (rows) × dimensions (cols) DataFrame.

        Reads scores stored in ``embed.varm`` by :meth:`calculate_interpretability_scores` and the
        per-dimension annotations from :meth:`set_latent_dimension_stats`, returning a tidy,
        ordered and labeled table for inspection.

        Parameters
        ----------
        embed
            AnnData of the latent space with scores in ``varm`` and the
            :meth:`set_latent_dimension_stats` annotations in ``var``.
        adata
            Data AnnData, used to take the gene names for the DataFrame columns.
        key
            ``varm`` score key to read (default ``"OOD_combined"``). With ``directional``, the
            ``_positive`` / ``_negative`` suffixes are appended automatically.
        directional
            If ``True`` (default), include both the positive and negative direction of each factor
            as separate columns (titles suffixed with ``+`` / ``-``).
        gene_symbols
            ``adata.var`` column to use for gene names; defaults to ``adata.var_names``.
        order_col
            ``embed.var`` column used to order the dimension columns (default ``"order"``).
        title_col
            ``embed.var`` column used to label the dimension columns (default ``"title"``).
        hide_vanished
            If ``True`` (default), drop factors (or directions) flagged as vanished by
            :meth:`set_latent_dimension_stats`.

        Returns
        -------
        A :class:`~pandas.DataFrame` of genes (rows) × dimensions (columns).
        """
        if directional:
            effect_data = np.concatenate(
                [embed.varm[key + "_positive"], embed.varm[key + "_negative"]]
            )
            var_info = (
                pd.concat([embed.var.assign(direction="+"), embed.var.assign(direction="-")])
                .assign(title=lambda df: df[title_col] + df["direction"])
                .reset_index(drop=True)
            )
            var_info["keep"] = (
                ~np.where(
                    var_info["direction"] == "+",
                    var_info["vanished_positive_direction"],
                    var_info["vanished_negative_direction"],
                )
                if hide_vanished
                else True
            )
        else:
            effect_data = embed.varm[key]
            var_info = embed.var.assign(title=lambda df: df[title_col]).assign(direction="")
            var_info["keep"] = ~var_info["vanished"] if hide_vanished else True

        gene_names = adata.var_names if gene_symbols is None else adata.var[gene_symbols]
        return (
            pd.DataFrame(effect_data, columns=gene_names, index=var_info["title"])
            .loc[var_info.query("keep == True").sort_values([order_col, "direction"])["title"]]
            .T
        )
