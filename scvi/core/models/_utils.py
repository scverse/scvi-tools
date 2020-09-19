import numpy as np
import pandas as pd
from scvi.core.utils import DifferentialComputation
from scvi._utils import track


def _de_core(
    adata,
    model_fn,
    groupby,
    group1,
    group2,
    idx1,
    idx2,
    all_stats,
    all_stats_fn,
    col_names,
    mode,
    batchid1,
    batchid2,
    delta,
    batch_correction,
    fdr,
    **kwargs
):
    """Internal function for DE interface."""
    if group1 is None and idx1 is None:
        group1 = adata.obs[groupby].cat.categories.tolist()

    if isinstance(group1, str):
        group1 = [group1]

    # make a temp obs key using indices
    temp_key = None
    if idx1 is not None:
        idx1 = np.asarray(idx1).ravel()
        g1_key = "one"
        obs_col = np.array(["None"] * adata.shape[0], dtype=str)
        obs_col[idx1] = g1_key
        group2 = None if idx2 is None else "two"
        if idx2 is not None:
            idx2 = np.asarray(idx2).ravel()
            obs_col[idx2] = group2
        temp_key = "_scvi_temp_de"
        adata.obs[temp_key] = obs_col
        groupby = temp_key
        group1 = [g1_key]

    df_results = []
    dc = DifferentialComputation(model_fn, adata)
    for g1 in track(
        group1,
        description="DE...",
    ):
        cell_idx1 = (adata.obs[groupby] == g1).ravel()
        if group2 is None:
            cell_idx2 = ~cell_idx1
        else:
            cell_idx2 = adata.obs[groupby] == group2

        all_info = dc.get_bayes_factors(
            cell_idx1,
            cell_idx2,
            mode=mode,
            delta=delta,
            batchid1=batchid1,
            batchid2=batchid2,
            use_observed_batches=not batch_correction,
            **kwargs,
        )

        if all_stats is True:
            genes_properties_dict = all_stats_fn(adata, cell_idx1, cell_idx2)
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=col_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        if mode == "change":
            res["is_de_fdr_{}".format(fdr)] = _fdr_de_prediction(
                res["proba_de"], fdr=fdr
            )
        if idx1 is None:
            g2 = "Rest" if group2 is None else group2
            res["comparison"] = "{} vs {}".format(g1, g2)
        df_results.append(res)

    if temp_key is not None:
        del adata.obs[temp_key]

    result = pd.concat(df_results, axis=0)

    return result


def _fdr_de_prediction(posterior_probas: np.ndarray, fdr: float = 0.05):
    """Compute posterior expected FDR and tag features as DE."""
    if not posterior_probas.ndim == 1:
        raise ValueError("posterior_probas should be 1-dimensional")
    sorted_genes = np.argsort(-posterior_probas)
    sorted_pgs = posterior_probas[sorted_genes]
    cumulative_fdr = (1.0 - sorted_pgs).cumsum() / (1.0 + np.arange(len(sorted_pgs)))
    d = (cumulative_fdr <= fdr).sum()
    pred_de_genes = sorted_genes[:d]
    is_pred_de = np.zeros_like(cumulative_fdr).astype(bool)
    is_pred_de[pred_de_genes] = True
    return is_pred_de
