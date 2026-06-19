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

        The effect formula depends on ``module.split_aggregation`` (``"logsumexp"`` or ``"mean"``).
        """
        for inference_outputs, generative_outputs in self.iterate_on_ae_output(
            adata=adata,
            dataloader=dataloader,
            deterministic=deterministic,
            **kwargs,
        ):
            # posterior mean, n_samples x n_latent
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
        """Effect of each split on reconstruction.

        Shape is ``(n_splits,)`` (or ``(2, n_splits)`` if directional) when aggregating over cells,
        else ``(n_cells, n_splits)`` (or ``(n_cells, 2, n_splits)``).
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

        Adds columns ``reconstruction_effect``, ``order``, ``max_value``, ``mean``, ``min``,
        ``max``, ``std``, ``std_abs``, ``title`` and ``vanished`` (+ directional vanished flags).
        """
        if embed.n_vars != self.module.n_split_latent:
            raise ValueError(
                "Per-dimension interpretability assumes one split per latent dimension "
                f"(n_split_latent == n_latent). Got n_split_latent={self.module.n_split_latent} "
                f"and {embed.n_vars} embedding dimensions."
            )
        if "original_dim_id" not in embed.var:
            embed.var["original_dim_id"] = np.arange(embed.var.shape[0])

        embed.var["reconstruction_effect"] = 0.0
        embed.var.loc[embed.var.sort_values("original_dim_id").index, "reconstruction_effect"] = (
            self.get_reconstruction_effect_of_each_split(adata=adata, dataloader=dataloader)
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

        Returns one ``(n_splits, n_genes)`` (or ``(2, n_splits, n_genes)`` if directional) array
        per requested aggregation (``max``, ``linear_weighted_mean``, ``exp_weighted_mean``).
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
                        weights = latent.unsqueeze(-1)  # n_samples x n_latent x 1
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
    ) -> dict[str, np.ndarray]:
        """Out-of-distribution per-split, per-gene effect scores via latent traversal.

        Traverses each latent dimension between its observed min/max over random batch/covariate
        combinations and returns ``min_possible``, ``max_possible`` and ``combined``
        log-fold-change style effects (``(n_splits, n_genes)`` or ``(2, n_splits, n_genes)`` if
        directional). Requires ``embed.var`` to contain ``min`` and ``max``
        (see :meth:`set_latent_dimension_stats`).
        """
        assert n_steps % 2 == 0, "n_steps must be even"
        dim_mins = np.minimum(embed.var["min"].values, 0.0)
        dim_maxs = np.maximum(embed.var["max"].values, 0.0)

        n_cat_total = [self.summary_stats.n_batch] + self._get_n_cats_per_cov()
        all_cat_combinations = np.asarray(
            list(itertools.product(*[range(n) for n in n_cat_total]))
        )
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

            min_per_split = effect_tensors.amin(dim=[0, 1])  # n_splits x n_genes
            max_per_split = effect_tensors.amax(dim=[0, 1])  # n_splits x n_genes
            directional_max_per_split = (
                effect_tensors.amax(dim=1) if directional else max_per_split
            )

            log_add_to_counts = torch.logsumexp(max_per_split, dim=[0, 1]) + np.log(
                add_to_counts / 1e6
            )
            add_min = log_add_to_counts.reshape(1, 1).expand(1, min_per_split.shape[1])
            add_max = log_add_to_counts.reshape(1, 1).expand(1, max_per_split.shape[1])
            lse_min = torch.logsumexp(
                torch.cat([min_per_split, add_min], dim=0), dim=0, keepdim=True
            )
            lse_max = torch.logsumexp(
                torch.cat([max_per_split, add_max], dim=0), dim=0, keepdim=True
            )
            max_possible = (
                torch.log(
                    torch.exp(lse_min)
                    - torch.exp(min_per_split)
                    + torch.exp(directional_max_per_split)
                )
                - lse_min
            )
            min_possible = torch.log(
                torch.exp(lse_max)
                - torch.exp(max_per_split)
                + torch.exp(directional_max_per_split)
            ) - torch.log(torch.exp(lse_max) - torch.exp(max_per_split) + torch.exp(min_per_split))
            combined = max_possible * min_possible

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
        """Compute in- and/or out-of-distribution interpretability scores.

        ``methods`` selects ``"IND"``, ``"OOD"`` (default) or ``"ALL"``. With ``inplace=True`` the
        scores are stored in ``embed.varm`` (keys
        ``"{method}_{aggregation}[_positive/_negative]"``); otherwise returned as a dict.
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

        Reads scores stored in ``embed.varm`` by :meth:`calculate_interpretability_scores`.
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
