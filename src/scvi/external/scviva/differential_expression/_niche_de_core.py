from collections.abc import Iterable as IterableClass

import pandas as pd
from rich import print

from scvi import REGISTRY_KEYS
from scvi.external.scviva import SCVIVA_REGISTRY_KEYS
from scvi.model.base._de_core import _fdr_de_prediction, _prepare_obs
from scvi.model.base._differential import DifferentialComputation
from scvi.utils import track

from ._de_utils import _get_nonzero_indices_from_rows, adjusted_nearest_neighbors
from ._marker_classifier import _gaussian_process_classifier
from ._results_dataclass import DifferentialExpressionResults


def _niche_de_core(
    adata_manager,
    model_fn,
    representation_fn,
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
    silent,
    ###### scVIVA specific ######
    radius=None,
    k_nn=None,
    lfc_select: str = "lfc_median",
    n_restarts_optimizer_gpc: int = 10,
    return_neighbors_idx: bool = True,
    **kwargs,
) -> DifferentialExpressionResults:
    """Internal function for DE interface."""
    adata = adata_manager.adata
    if group1 is None and idx1 is None:
        group1 = adata.obs[groupby].astype("category").cat.categories.tolist()
        if len(group1) == 1:
            raise ValueError("Only a single group in the data. Can't run DE on a single group.")

    if not isinstance(group1, IterableClass) or isinstance(group1, str):
        group1 = [group1]

    # make a temp obs key using indices
    temp_key = None
    if idx1 is not None:
        obs_col, group1, group2 = _prepare_obs(idx1, idx2, adata)
        temp_key = "_scvi_temp_de"
        adata.obs[temp_key] = obs_col
        groupby = temp_key

    cell_samples = adata_manager.get_from_registry(SCVIVA_REGISTRY_KEYS.SAMPLE_KEY)
    cell_labels = adata_manager.get_from_registry(REGISTRY_KEYS.LABELS_KEY)
    cell_coordinates = adata_manager.get_from_registry(SCVIVA_REGISTRY_KEYS.CELL_COORDINATES_KEY)

    # don't compute adjusted nearest neighbors if already computed
    if "adjusted_A" in adata.uns.keys():
        A = adata.uns["adjusted_A"]
    else:
        print("Computing adjusted nearest neighbors...")
        A = adjusted_nearest_neighbors(
            cell_samples=cell_samples,
            cell_coordinates=cell_coordinates,
            cell_labels=cell_labels,
            radius=radius,
            k_nn=k_nn,
            return_sparse=True,
        )

        adata.uns["adjusted_A"] = A

    print("Computing DE...")
    if group2 is not None:
        comparisons = [
            "group1_group2",
            "group1_neighbors1",
            "neighbors1_group2",
            "neighbors1_neighbors2",
        ]
        DE_results = {comparison: [] for comparison in comparisons}

    else:
        comparisons = ["group1_group2", "group1_neighbors1", "neighbors1_group2"]
        DE_results = {comparison: [] for comparison in comparisons}
    dc = DifferentialComputation(model_fn, representation_fn, adata_manager)
    for g1 in track(
        group1,
        description="DE...",
        disable=silent,
    ):
        cell_idx1 = (adata.obs[groupby] == g1).to_numpy().ravel()
        neighbors_idx1 = _get_nonzero_indices_from_rows(A, cell_idx1)

        if group2 is None:
            cell_idx2 = ~cell_idx1

            DE_indices = {
                "group1_group2": [cell_idx1, cell_idx2],
                "group1_neighbors1": [cell_idx1, neighbors_idx1],
                "neighbors1_group2": [neighbors_idx1, cell_idx2],
            }
            DE_group_names = {
                "group1_group2": [g1, "Rest"],
                "group1_neighbors1": [g1, f"{g1}_neighbors"],
                "neighbors1_group2": [f"{g1}_neighbors", "Rest"],
            }

        else:
            cell_idx2 = (adata.obs[groupby] == group2).to_numpy().ravel()
            neighbors_idx2 = _get_nonzero_indices_from_rows(A, cell_idx2)
            DE_indices = {
                "group1_group2": [cell_idx1, cell_idx2],
                "group1_neighbors1": [cell_idx1, neighbors_idx1],
                "neighbors1_group2": [neighbors_idx1, cell_idx2],
                "neighbors1_neighbors2": [neighbors_idx1, neighbors_idx2],
            }
            DE_group_names = {
                "group1_group2": [g1, group2],
                "group1_neighbors1": [g1, f"{g1}_neighbors"],
                "neighbors1_group2": [f"{g1}_neighbors", group2],
                "neighbors1_neighbors2": [f"{g1}_neighbors", f"{group2}_neighbors"],
            }

        # Ensure fdr and delta are lists of the correct length
        if isinstance(fdr, list):
            assert len(fdr) == len(comparisons), (
                f"Mismatch: len(fdr)={len(fdr)}, expected={len(comparisons)}"
            )
        else:
            fdr = [fdr] * len(comparisons)  # Convert to list of same value

        if isinstance(delta, list):
            assert len(delta) == len(comparisons), (
                f"Mismatch: len(delta)={len(delta)}, expected={len(comparisons)}"
            )
        else:
            delta = [delta] * len(comparisons)  # Convert to list of same value

        # Assign FDR and Delta values
        DE_group_fdr = dict(zip(comparisons, fdr, strict=False))
        DE_group_delta = dict(zip(comparisons, delta, strict=False))

        for comparison, [cell_idx_a, cell_idx_b] in DE_indices.items():
            print(f"Running DE for {comparison}")

            all_info = dc.get_bayes_factors(
                cell_idx_a,
                cell_idx_b,
                mode=mode,
                delta=DE_group_delta[comparison],
                batchid1=batchid1,
                batchid2=batchid2,
                use_observed_batches=not batch_correction,
                test_mode="three",
                **kwargs,
            )

            if all_stats is True:
                genes_properties_dict = all_stats_fn(adata_manager, cell_idx_a, cell_idx_b)
                all_info = {**all_info, **genes_properties_dict}

            res = pd.DataFrame(all_info, index=col_names)
            sort_key = "proba_de" if mode == "change" else "bayes_factor"
            res = res.sort_values(by=sort_key, ascending=False)
            if mode == "change":
                res[f"is_de_fdr_{DE_group_fdr[comparison]}"] = _fdr_de_prediction(
                    res["proba_de"], fdr=DE_group_fdr[comparison]
                )
            if idx1 is None:
                g1_name = DE_group_names[comparison][0]
                g2_name = DE_group_names[comparison][1]
                res["comparison"] = f"{g1_name} vs {g2_name}"
                res["group1"] = g1_name
                res["group2"] = g2_name
            DE_results[comparison].append(res)

    if temp_key is not None:
        del adata.obs[temp_key]

    DE_results["group1_group2"] = pd.concat(DE_results["group1_group2"], axis=0)
    idx_g1_g2 = DE_results["group1_group2"].index

    if group2 is None:
        group2 = "Rest"

    for groups in list(DE_results.keys())[1:]:
        group_DE_result = DE_results[groups]
        DE_results[groups] = pd.concat(group_DE_result, axis=0).reindex(idx_g1_g2)

    # fit the classifier
    lfc_g1_g2 = DE_results["group1_group2"][lfc_select]
    lfc_n1_g2 = DE_results["neighbors1_group2"][lfc_select]
    fdr_g1_n1 = DE_results["group1_neighbors1"][
        f"is_de_fdr_{DE_group_fdr['group1_neighbors1']}"
    ].copy()
    fdr_g1_n1.loc[DE_results["group1_neighbors1"][lfc_select] < 0] = False

    if fdr_g1_n1.sum() == 0:
        raise ValueError("No DE genes found between group1 and neighbors1.")

    print("Computing g1 confidence scores...")
    gpc_ = _gaussian_process_classifier(
        lfc_g1_g2,
        lfc_n1_g2,
        fdr_g1_n1,
        n_restarts_optimizer=n_restarts_optimizer_gpc,
        restrict_to_upregulated=True,
    )

    for groups in list(DE_results.keys()):
        DE_results[groups]["proba_de_g1_n1"] = gpc_.gene_probas_

    return DifferentialExpressionResults(
        gpc=gpc_,
        g1_g2=DE_results["group1_group2"],
        g1_n1=DE_results["group1_neighbors1"],
        n1_g2=DE_results["neighbors1_group2"],
        n1_n2=DE_results["neighbors1_neighbors2"]
        if "neighbors1_neighbors2" in DE_results
        else None,
        n1_index=neighbors_idx1 if return_neighbors_idx else None,
        n2_index=neighbors_idx2 if return_neighbors_idx and group2 != "Rest" else None,
    )
