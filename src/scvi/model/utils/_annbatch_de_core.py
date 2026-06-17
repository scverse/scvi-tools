import inspect
from collections.abc import Callable
from collections.abc import Iterable as IterableClass

import numpy as np
import pandas as pd
import torch

from scvi.model.base._differential import (
    describe_continuous_distrib,
    estimate_delta,
    estimate_pseudocounts_offset,
    pairs_sampler,
)
from scvi.utils import track


def _as_numpy_array(x):
    """Convert dataloader tensors and sparse matrices to dense numpy arrays."""
    if isinstance(x, torch.Tensor):
        if x.is_sparse:
            x = x.to_dense()
        return x.detach().cpu().numpy()
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def _stream_raw_counts_properties_from_dataloader(
    dataloader,
    mask1: np.ndarray,
    mask2: np.ndarray,
    n_vars: int,
    var_idx: list[int] | np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Compute raw count statistics for streamed annbatch data."""
    if var_idx is not None:
        var_idx = np.asarray(var_idx)
        n_vars = len(var_idx)

    batch_cursor = 0
    sum1 = np.zeros(n_vars)
    sum2 = np.zeros(n_vars)
    nonz1 = np.zeros(n_vars)
    nonz2 = np.zeros(n_vars)
    norm_sum1 = np.zeros(n_vars)
    norm_sum2 = np.zeros(n_vars)
    n1 = 0
    n2 = 0

    for batch in dataloader:
        x = _as_numpy_array(batch["X"])
        batch_size = x.shape[0]
        sl = slice(batch_cursor, batch_cursor + batch_size)
        local1 = np.asarray(mask1[sl])
        local2 = np.asarray(mask2[sl])
        if var_idx is not None:
            x = x[:, var_idx]
        if local1.any():
            x1 = x[local1]
            sum1 += x1.sum(axis=0)
            nonz1 += (x1 != 0).mean(axis=0) * x1.shape[0]
            scaling = 1 / np.asarray(x1.sum(axis=1)).ravel()
            scaling *= 1e4
            norm_sum1 += (x1 * scaling[:, None]).sum(axis=0)
            n1 += x1.shape[0]
        if local2.any():
            x2 = x[local2]
            sum2 += x2.sum(axis=0)
            nonz2 += (x2 != 0).mean(axis=0) * x2.shape[0]
            scaling = 1 / np.asarray(x2.sum(axis=1)).ravel()
            scaling *= 1e4
            norm_sum2 += (x2 * scaling[:, None]).sum(axis=0)
            n2 += x2.shape[0]
        batch_cursor += batch_size

    return {
        "raw_mean1": sum1 / max(n1, 1),
        "raw_mean2": sum2 / max(n2, 1),
        "non_zeros_proportion1": nonz1 / max(n1, 1),
        "non_zeros_proportion2": nonz2 / max(n2, 1),
        "raw_normalized_mean1": norm_sum1 / max(n1, 1),
        "raw_normalized_mean2": norm_sum2 / max(n2, 1),
    }


def _sample_normalized_expression_from_dataloader(
    model_fn: Callable,
    mask1: np.ndarray,
    mask2: np.ndarray,
    n_samples_overall: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample normalized expression from an annbatch-backed model function."""
    target_samples = 5000 if n_samples_overall is None else int(n_samples_overall)
    min_cells = min(int(np.sum(mask1)), int(np.sum(mask2)))
    if min_cells == 0:
        raise ValueError("Both groups must contain at least one cell.")
    n_samples_per_cell = target_samples // min_cells + 1

    expr = model_fn(n_samples=n_samples_per_cell)
    expr = np.asarray(expr)
    if expr.ndim == 2:
        expr = expr[None, ...]

    def _select_samples(mask):
        selected = expr[:, mask, :].reshape(-1, expr.shape[-1])
        sample_idx = np.random.choice(
            selected.shape[0],
            target_samples,
            replace=True,
        )
        return selected[sample_idx]

    return _select_samples(mask1), _select_samples(mask2)


def _annbatch_bayes_factors_from_samples(
    scales_1: np.ndarray,
    scales_2: np.ndarray,
    mean1: np.ndarray,
    mean2: np.ndarray,
    mode: str,
    delta: float | None,
    raw_stats: dict[str, np.ndarray] | None,
    use_permutation: bool,
    m_permutation: int,
    change_fn: str | Callable | None,
    m1_domain_fn: Callable | None,
    pseudocounts: float | None,
    threshold_counts: float,
    test_mode: str,
    cred_interval_lvls: list[float] | np.ndarray | None,
) -> dict[str, np.ndarray]:
    """Compute DE Bayes-factor fields from annbatch posterior samples."""
    eps = 1e-8

    if mode == "vanilla":
        scales_1, scales_2 = pairs_sampler(
            scales_1,
            scales_2,
            use_permutation=use_permutation,
            m_permutation=m_permutation,
        )
        proba_m1 = np.mean(scales_1 > scales_2, axis=0)
        proba_m2 = 1.0 - proba_m1
        return {
            "proba_m1": proba_m1,
            "proba_m2": proba_m2,
            "bayes_factor": np.log(proba_m1 + eps) - np.log(proba_m2 + eps),
            "scale1": mean1,
            "scale2": mean2,
        }

    if mode != "change":
        raise NotImplementedError(f"Mode {mode} not recognized")

    scales_1, scales_2 = pairs_sampler(
        scales_1,
        scales_2,
        use_permutation=use_permutation,
        m_permutation=m_permutation,
    )
    if pseudocounts is None:
        if raw_stats is None:
            raise ValueError("Raw statistics are required to estimate pseudocounts.")
        pseudocounts = estimate_pseudocounts_offset(
            scales_a=scales_1,
            scales_b=scales_2,
            where_zero_a=raw_stats["raw_mean1"] < threshold_counts,
            where_zero_b=raw_stats["raw_mean2"] < threshold_counts,
            quantile=0.9,
        )

    def lfc(x, y, pseudocounts=pseudocounts):
        return np.log2(x + pseudocounts) - np.log2(y + pseudocounts)

    change_fn_ = lfc if change_fn == "log-fold" or change_fn is None else change_fn
    if not callable(change_fn_):
        raise ValueError("'change_fn' attribute not understood")

    m1_domain_fn_ = m1_domain_fn
    if m1_domain_fn_ is None:

        def m1_domain_fn_(samples):
            delta_ = delta if delta is not None else estimate_delta(lfc_means=samples.mean(axis=0))
            samples_plus = samples >= delta_
            samples_minus = samples < -delta_
            return samples_plus, samples_minus

    change_fn_specs = inspect.getfullargspec(change_fn_)
    domain_fn_specs = inspect.getfullargspec(m1_domain_fn_)
    if (len(change_fn_specs.args) != 3) or (len(domain_fn_specs.args) != 1):
        raise ValueError(
            "change_fn should take exactly three parameters as inputs; m1_domain_fn one parameter."
        )
    try:
        change_distribution = change_fn_(scales_1, scales_2, pseudocounts)
        is_de_plus, is_de_minus = m1_domain_fn_(change_distribution)
        delta_ = (
            estimate_delta(lfc_means=change_distribution.mean(axis=0)) if delta is None else delta
        )
    except TypeError as err:
        raise TypeError(
            "change_fn or m1_domain_fn have has wrong properties."
            "Please ensure that these functions have the right signatures and"
            "outputs and that they can process numpy arrays"
        ) from err

    proba_m1 = np.mean(is_de_plus, axis=0)
    proba_m2 = np.mean(is_de_minus, axis=0)
    proba_de = proba_m1 + proba_m2 if test_mode == "two" else np.maximum(proba_m1, proba_m2)
    change_distribution_props = describe_continuous_distrib(
        samples=change_fn_(scales_1, scales_2, 1e-3 * pseudocounts),
        credible_intervals_levels=cred_interval_lvls,
    )
    change_distribution_props = {
        "lfc_" + key: val for (key, val) in change_distribution_props.items()
    }
    return {
        "proba_de": proba_de,
        "proba_not_de": 1.0 - proba_de,
        "bayes_factor": np.log(proba_de + eps) - np.log(1.0 - proba_de + eps),
        "scale1": mean1,
        "scale2": mean2,
        "pseudocounts": pseudocounts,
        "delta": delta_,
        **change_distribution_props,
    }


def _de_core_for_annbatch(
    dataloader,
    model_fn: Callable,
    obs: pd.DataFrame,
    groupby,
    group1,
    group2,
    idx1,
    idx2,
    all_stats,
    col_names,
    mode,
    delta,
    fdr,
    silent,
    **kwargs,
):
    """Internal DE interface for annbatch-backed models."""
    if groupby is None and idx1 is None:
        raise ValueError("Must provide `groupby` or `idx1` when using a dataloader.")

    col_names = np.asarray(col_names)
    n_samples_overall = kwargs.get("n_samples_overall", 5000)
    use_permutation = kwargs.get("use_permutation", False)
    m_permutation = kwargs.get("m_permutation", 10000)
    change_fn = kwargs.get("change_fn", None)
    m1_domain_fn = kwargs.get("m1_domain_fn", None)
    pseudocounts = kwargs.get("pseudocounts", None)
    threshold_counts = kwargs.get("threshold_counts", 0.01)
    test_mode = kwargs.get("test_mode", "three")
    cred_interval_lvls = kwargs.get("cred_interval_lvls", None)

    def _to_indices(mask):
        if mask is None:
            return None
        mask = np.asarray(mask)
        if mask.dtype == bool:
            return np.where(mask)[0]
        return mask

    if idx1 is not None:
        idx1 = _to_indices(idx1)
        idx2 = _to_indices(idx2)
    elif group1 is None:
        group1 = obs[groupby].astype("category").cat.categories.tolist()
        if len(group1) == 1:
            raise ValueError("Only a single group in the data. Can't run DE on a single group.")

    if group1 is not None and (not isinstance(group1, IterableClass) or isinstance(group1, str)):
        group1 = [group1]

    df_results = []
    groups = [None] if idx1 is not None else group1
    for g1 in track(groups, description="DE...", disable=silent):
        if idx1 is None:
            cell_idx1 = (obs[groupby] == g1).to_numpy().ravel()
            cell_idx2 = (
                ~cell_idx1 if group2 is None else (obs[groupby] == group2).to_numpy().ravel()
            )
        else:
            cell_idx1 = np.zeros(len(obs), dtype=bool)
            cell_idx1[np.asarray(idx1)] = True
            if idx2 is None:
                cell_idx2 = ~cell_idx1
            else:
                cell_idx2 = np.zeros(len(obs), dtype=bool)
                cell_idx2[np.asarray(idx2)] = True

        scales_1, scales_2 = _sample_normalized_expression_from_dataloader(
            model_fn,
            cell_idx1,
            cell_idx2,
            n_samples_overall,
        )
        mean1 = scales_1.mean(axis=0)
        mean2 = scales_2.mean(axis=0)
        raw_stats = (
            _stream_raw_counts_properties_from_dataloader(
                dataloader,
                cell_idx1,
                cell_idx2,
                len(col_names),
            )
            if all_stats or (mode == "change" and pseudocounts is None)
            else None
        )
        all_info = _annbatch_bayes_factors_from_samples(
            scales_1=scales_1,
            scales_2=scales_2,
            mean1=mean1,
            mean2=mean2,
            mode=mode,
            delta=delta,
            raw_stats=raw_stats,
            use_permutation=use_permutation,
            m_permutation=m_permutation,
            change_fn=change_fn,
            m1_domain_fn=m1_domain_fn,
            pseudocounts=pseudocounts,
            threshold_counts=threshold_counts,
            test_mode=test_mode,
            cred_interval_lvls=cred_interval_lvls,
        )

        if all_stats:
            all_info = {**all_info, **raw_stats}
        res = pd.DataFrame(all_info, index=col_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        if mode == "change":
            from scvi.model.base._de_core import _fdr_de_prediction

            res[f"is_de_fdr_{fdr}"] = _fdr_de_prediction(res["proba_de"], fdr=fdr)
        if idx1 is None:
            g2 = "Rest" if group2 is None else group2
            res["comparison"] = f"{g1} vs {g2}"
            res["group1"] = g1
            res["group2"] = g2
        df_results.append(res)

    return pd.concat(df_results, axis=0)
