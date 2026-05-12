from functools import partial
from typing import Literal, Optional, Union

import anndata
import numpy as np
import pandas as pd

from visionpy import data_accessor, hotspot_acc, vision_acc

from ..preprocessing.anndata import setup_anndata
from .cell_communication import (
    compute_cell_communication_np_fast,
    compute_cell_communication_p_fast,
    compute_gene_pairs,
)
from .knn import (
    compute_neighbors,
    compute_neighbors_from_distances,
    compute_weights,
    filter_distances,
    get_indexes,
    object_scalar,
)
from .local_autocorrelation import compute_local_autocorrelation, compute_local_autocorrelation_fast
from .local_correlation import compute_local_correlation_fast
from .modules import calculate_module_scores, compute_modules
from .signature import compute_signatures_anndata


def start_vision(
    adata: Union[str, anndata.AnnData],
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    distances_obsp_key: Optional[str] = None,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    weighted_graph: Optional[bool] = False,
    neighborhood_radius: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    neighborhood_factor: Optional[int] = 3,
    sample_key: Optional[str] = None,
    obs_df_scores: Optional[bool] = False,
    one_vs_all_obs_cols: Optional[bool] = False,
    one_vs_all_signatures: Optional[bool] = False,
    gene_score_per_signature: Optional[bool] = False,
    scores_only: Optional[bool] = False,
):
    """Start VISION.

    Parameters
    ----------
    adata
        AnnData object.
    norm_data_key
        Key for layer with log library size normalized data. If `None` (default), uses `adata.X`.
    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use neighbors stored in `adata`. If no neighbors have been previously computed an error will be raised.
    distances_obsp_key
        Distances encoding cell-cell similarities directly. Shape is (cells x cells). Input is key in `adata.obsp`.
    signature_varm_key
        Key in `adata.varm` for signatures. If `None` (default), no signatures. Matrix should encode positive genes with 1, negative genes with -1, and all other genes with 0.
    signature_names_uns_key
        Key in `adata.uns` for signature names. If `None`, attempts to read columns if `signature_varm_key` is a pandas DataFrame. Otherwise, uses `Signature_1`, `Signature_2`, etc.
    weighted_graph
        Whether or not to create a weighted graph.
    neighborhood_radius
        Neighborhood radius.
    n_neighbors
        Neighborhood size.
    neighborhood_factor
        Used when creating a weighted graph.  Sets how quickly weights decay relative to the distances within the neighborhood. The weight for a cell with a distance d will decay as exp(-d^2/D) where D is the distance to the `n_neighbors`/`neighborhood_factor`-th neighbor.
    sample_key
        Sample information in case the data contains different samples or samples from different conditions. Input is key in `adata.obs`.
    obs_df_scores
        Boolean variable indicating whether to compute observation scores or not.
    one_vs_all_obs_cols
        Boolean variable indicating whether to compute one vs all DE analysis of the numerical variables for every categorical variable or not.
    one_vs_all_signatures
        Boolean variable indicating whether to compute one vs all DE analysis of the signature scores for every categorical variable or not.
    gene_score_per_signature
        Boolean variable indicating whether to compute the correlation between gene expression and signature scores for every gene-signature pair or not.

    """
    if isinstance(adata, str):
        adata = anndata.read(str)

    if scores_only is False:

        # compute neighbors and weights
        if compute_neighbors_on_key is not None:
            print("Computing the neighborhood graph...")
            compute_neighbors(
                adata=adata,
                compute_neighbors_on_key=compute_neighbors_on_key,
                n_neighbors=n_neighbors,
                neighborhood_radius=neighborhood_radius,
                sample_key=sample_key,
            )
        else:
            if distances_obsp_key is not None and distances_obsp_key in adata.obsp:
                print("Computing the neighborhood graph from distances...")
                compute_neighbors_from_distances(
                    adata,
                    distances_obsp_key,
                    n_neighbors,
                    sample_key,
                )

        if 'weights' not in adata.obsp and 'distances' in adata.obsp:
            print("Computing the weights...")
            compute_weights(
                adata,
                weighted_graph,
                neighborhood_factor,
            )

    vision_acc.adata = adata
    vision_acc.norm_data_key = norm_data_key
    vision_acc.compute_neighbors_on_key = compute_neighbors_on_key
    vision_acc.distances_obsp_key = distances_obsp_key
    vision_acc.signature_varm_key = signature_varm_key
    vision_acc.signature_names_uns_key = signature_names_uns_key
    vision_acc.weighted_graph = weighted_graph
    vision_acc.neighborhood_radius = neighborhood_radius
    vision_acc.n_neighbors = n_neighbors
    vision_acc.neighborhood_factor = neighborhood_factor
    vision_acc.sample_key = sample_key

    if obs_df_scores is True:
        vision_acc.compute_obs_df_scores()

    if one_vs_all_obs_cols is True:
        vision_acc.compute_one_vs_all_obs_cols()

    # compute signatures
    if signature_varm_key is not None:
        adata.obsm["vision_signatures"] = compute_signatures_anndata(
            adata,
            norm_data_key,
            signature_varm_key,
            signature_names_uns_key,
        )

        if scores_only is False:
            vision_acc.compute_signature_scores()

        if one_vs_all_signatures is True:
            vision_acc.compute_one_vs_all_signatures()

    if gene_score_per_signature is True:
        vision_acc.compute_gene_score_per_signature()


def start_hotspot(
    adata: Union[str, anndata.AnnData],
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    layer_key: Optional[Union[Literal["use_raw"], str]] = None,
    model: Optional[str] = None,
    compute_neighbors_on_key: Optional[str] = None,
    distances_obsp_key: Optional[str] = None,
    signature_varm_key: Optional[str] = None,
    weighted_graph: Optional[bool] = True,
    neighborhood_radius: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    neighborhood_factor: Optional[int] = 3,
    sample_key: Optional[str] = None,
    min_gene_threshold: Optional[int] = 15,
    core_only: Optional[bool] = False,
    fdr_threshold: Optional[float] = 0.05,
    jobs: Optional[int] = None,
):
    """Start Hotspot.

    Parameters
    ----------
    adata
        AnnData object.
    norm_data_key
        Key for layer with log library size normalized data. If `None` (default), uses `adata.X`.
    layer_key
        Key(s) in adata.layers with count data. If `None` (default), uses `adata.X`.
    model
        Specifies the null model to use for gene expression.
        Valid choices are:

            - 'danb': Depth-Adjusted Negative Binomial
            - 'bernoulli': Models probability of detection
            - 'normal': Depth-Adjusted Normal
            - 'none': Assumes data has been pre-standardized

    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use neighbors stored in `adata`. If no neighbors have been previously computed an error will be raised.
    distances_obsp_key
        Distances encoding cell-cell similarities directly. Shape is (cells x cells). Input is key in `adata.obsp`.
    weighted_graph
        Whether or not to create a weighted graph.
    neighborhood_radius
        Neighborhood radius.
    n_neighbors
        Neighborhood size.
    neighborhood_factor
        Used when creating a weighted graph.  Sets how quickly weights decay relative to the distances within the neighborhood. The weight for a cell with a distance d will decay as exp(-d^2/D) where D is the distance to the `n_neighbors`/`neighborhood_factor`-th neighbor.
    sample_key
        Sample information in case the data contains different samples or samples from different conditions. Input is key in `adata.obs`.
    jobs
        Number of parallel jobs to run.

    """
    if isinstance(adata, str):
        adata = anndata.read(str)

    # compute neighbors and weights
    if compute_neighbors_on_key is not None:
        print("Computing the neighborhood graph...")
        compute_neighbors(
            adata=adata,
            compute_neighbors_on_key=compute_neighbors_on_key,
            n_neighbors=n_neighbors,
            neighborhood_radius=neighborhood_radius,
            sample_key=sample_key,
        )
    else:
        if distances_obsp_key is not None and distances_obsp_key in adata.obsp:
            print("Computing the neighborhood graph from distances...")
            compute_neighbors_from_distances(
                adata,
                distances_obsp_key,
                n_neighbors,
                sample_key,
            )

    # distances = adata.obsp['distances']
    # distances = distances.tocsr()

    # dist = list(map(partial(filter_distances), distances))
    # dist = pd.concat(dist, axis=1).T
    # dist.index = adata.obs.index

    # ind = list(map(partial(get_indexes), distances))
    # ind = pd.concat(ind, axis=1).T
    # ind.index = adata.obs.index

    if 'weights' not in adata.obsp and 'distances' in adata.obsp:
        print("Computing the weights...")
        compute_weights(
            adata,
            weighted_graph,
            neighborhood_factor,
        )

    hotspot_acc.adata = adata
    hotspot_acc.norm_data_key = norm_data_key
    hotspot_acc.layer_key = layer_key
    hotspot_acc.model = model
    hotspot_acc.compute_neighbors_on_key = compute_neighbors_on_key
    hotspot_acc.distances_obsp_key = distances_obsp_key
    hotspot_acc.signature_varm_key = signature_varm_key
    hotspot_acc.weighted_graph = weighted_graph
    hotspot_acc.neighborhood_radius = neighborhood_radius
    hotspot_acc.n_neighbors = n_neighbors
    hotspot_acc.neighborhood_factor = neighborhood_factor
    hotspot_acc.sample_key = sample_key
    hotspot_acc.min_gene_threshold = min_gene_threshold
    hotspot_acc.core_only = core_only
    hotspot_acc.fdr_threshold = fdr_threshold
    hotspot_acc.jobs = jobs

    # compute gene local autocorrelation
    print("Computing local autocorrelation...")
    compute_local_autocorrelation_fast(
        adata=adata,
        layer_key=layer_key,
        model=model,
    )

    # Select the genes with significant autocorrelation
    gene_autocorrelation_results = adata.uns['gene_autocorrelation_results']
    genes = gene_autocorrelation_results.loc[gene_autocorrelation_results.FDR < 0.05].sort_values('Z', ascending=False).index

    # compute local correlation
    print("Computing local correlation...")
    compute_local_correlation_fast(
        adata,
        layer_key,
        model,
        genes,
    )

    # compute modules
    print("Computing modules...")
    compute_modules(
        adata,
        min_gene_threshold=min_gene_threshold,
        core_only=core_only,
        fdr_threshold=fdr_threshold,
    )

    # compute modules
    print("Computing module scores...")
    calculate_module_scores(
        adata,
        layer_key,
        model,
    )

    if ("vision_signatures" in adata.obsm) and (len(adata.uns["gene_modules_dict"].keys()) > 0):
        print("Computing signature-module enrichment...")
        hotspot_acc.compute_sig_mod_enrichment()

        adata.obsm["signature_modules_overlap"] = compute_signatures_anndata(
            adata,
            norm_data_key,
            signature_varm_key='signatures_overlap',
            signature_names_uns_key=None,
        )


def start_CCC_analysis(
    adata: Union[str, anndata.AnnData],
    layer_key: Optional[Union[Literal["use_raw"], str]] = None,
    model: Optional[str] = None,
    compute_neighbors_on_key: Optional[str] = None,
    distances_obsp_key: Optional[str] = None,
    deconv_data: Optional[bool] = False,
    cell_type_list: Optional[list] = None,
    cell_type_key: Optional[str] = None,
    cell_type_pairs: Optional[list] = None,
    database_varm_key: Optional[str] = None,
    weighted_graph: Optional[bool] = True,
    neighborhood_radius: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    neighborhood_factor: Optional[int] = 3,
    spot_diameter: Optional[int] = None,
    sample_key: Optional[str] = None,
    autocorrelation_filt: Optional[bool] = False,
    expression_filt: Optional[bool] = False,
    de_filt: Optional[bool] = False,
    test: Optional[str] = None,
    jobs: Optional[int] = None,
    run_CCC_pipeline: Optional[bool] = True,
    return_all_adatas: Optional[bool] = False,
):
    """Start the CCC analysis.

    Parameters
    ----------
    adata
        AnnData object.
    layer_key
        Key(s) in adata.layers with count data. If `None` (default), uses `adata.X`
    model
        Specifies the null model to use for gene expression.
        Valid choices are:

            - 'danb': Depth-Adjusted Negative Binomial
            - 'bernoulli': Models probability of detection
            - 'normal': Depth-Adjusted Normal
            - 'none': Assumes data has been pre-standardized

    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use neighbors stored in `adata`. If no neighbors have been previously computed an error will be raised.
    distances_obsp_key
        Distances encoding cell-cell similarities directly. Shape is (cells x cells). Input is key in `adata.obsp`.
    deconv_data
        Boolean variable indicating whether the dataset has been deconvolved or not.
    cell_type_list
        Cell type or cluster information for the cell-cell communication analysis. Input is a list of keys in `adata.layers`.
    cell_type_key
        Cell type or cluster information for the cell-cell communication analysis. Input is key in `adata.obs`.
    database_uns_key
        Gene pair database to filter interactions. If not provided, all gene pairs from the AnnData are considered. qw
    weighted_graph
        Whether or not to create a weighted graph.
    neighborhood_radius
        Neighborhood radius.
    n_neighbors
        Neighborhood size.
    neighborhood_factor
        Used when creating a weighted graph.  Sets how quickly weights decay relative to the distances within the neighborhood. The weight for a cell with a distance d will decay as exp(-d^2/D) where D is the distance to the `n_neighbors`/`neighborhood_factor`-th neighbor.
    spot_diameter
        Spot diameter of the spatial transcriptomics technology if deconvolution is performed.
    sample_key
        Sample information in case the data contains different samples or samples from different conditions. Input is key in `adata.obs`.
    autocorrelation_filt
        Boolean variable indicating whether to perform the autocorrelation filter or not.
    expression_filt
        Boolean variable indicating whether to perform the expression filter or not.
    de_filt
        Boolean variable indicating whether to perform the differential expression filter or not.
    test
        Whether to compute cell-cell communication using a parametric or a non-parametric test. Input should be either 'parametric' or 'non-parametric'.
    jobs
        Number of parallel jobs to run.
    run_CCC_pipeline
        Boolean variable indicating whether to run the cell-cell communication pipeline or not.

    """
    if isinstance(adata, str):
        adata = anndata.read(str)

    # setup AnnData
    if deconv_data is True:
        old_adata = adata.copy()
        adata = setup_anndata(
            adata,
            cell_type_list,
            compute_neighbors_on_key,
            cell_type_key,
            database_varm_key,
            sample_key,
            )
        layer_key=None

    # compute neighbors and weights
    if compute_neighbors_on_key is not None:
        print("Computing the neighborhood graph...")
        compute_neighbors(
            adata=adata,
            compute_neighbors_on_key=compute_neighbors_on_key,
            n_neighbors=n_neighbors,
            neighborhood_radius=neighborhood_radius,
            spot_diameter=spot_diameter,
            sample_key=sample_key,
            deconv_data=deconv_data,
        )
    else:
        if distances_obsp_key is not None and distances_obsp_key in adata.obsp:
            print("Computing the neighborhood graph from distances...")
            compute_neighbors_from_distances(
                adata=adata,
                distances_obsp_key=distances_obsp_key,
                spot_diameter=spot_diameter,
                sample_key=sample_key,
                deconv_data=deconv_data,
            )

    if 'weights' not in adata.obsp and 'distances' in adata.obsp:
        print("Computing the weights...")
        compute_weights(
            adata,
            weighted_graph,
            neighborhood_factor,
        )

    data_accessor.adata = adata
    data_accessor.layer_key = layer_key
    data_accessor.model = model
    data_accessor.compute_neighbors_on_key = compute_neighbors_on_key
    data_accessor.distances_obsp_key = distances_obsp_key
    data_accessor.deconv_data = deconv_data
    data_accessor.cell_type_list = cell_type_list
    data_accessor.cell_type_key = cell_type_key
    data_accessor.cell_type_pairs = cell_type_pairs
    data_accessor.database_varm_key = database_varm_key
    data_accessor.weighted_graph = weighted_graph
    data_accessor.n_neighbors = n_neighbors
    data_accessor.neighborhood_factor = neighborhood_factor
    data_accessor.neighborhood_radius = neighborhood_radius
    data_accessor.spot_diameter = spot_diameter
    data_accessor.sample_key = sample_key
    data_accessor.autocorrelation_filt = autocorrelation_filt
    data_accessor.expression_filt = expression_filt
    data_accessor.de_filt = de_filt
    data_accessor.test = test
    data_accessor.jobs = jobs

    if run_CCC_pipeline is True:
        print('Running the cell-cell communication pipeline...')

        # compute gene local autocorrelation
        if ('gene_autocorrelation_results' not in adata.uns) and (autocorrelation_filt is True):
            compute_local_autocorrelation_fast(
                adata,
                layer_key,
                database_varm_key,
                model,
                jobs,
            )

        if (expression_filt is True) or (de_filt is True):
            data_accessor.filter_genes()

        compute_gene_pairs(
            adata,
            database_varm_key,
            layer_key,
            cell_type_key,
            autocorrelation_filt,
            expression_filt,
            de_filt,
            cell_type_pairs
        )

        # compute cell-cell communication (assess if it would be better to call these functions from the jupyter notebook)
        if test == 'parametric':
            compute_cell_communication_p_fast(
                adata,
                layer_key,
                database_varm_key,
                model,
                cell_type_key,
                cell_type_list,
            )
        elif test == 'non-parametric':
            compute_cell_communication_np_fast(
                adata,
                layer_key,
                database_varm_key,
                cell_type_key,
            )
        else:
            raise Exception(f"Invalid cell-cell communication test: {test}. It should be either 'parametric' or 'non-parametric'")

        if deconv_data is True:
            old_adata.uns["cell_communication_df"] = adata.uns["cell_communication_df"]
            old_adata.uns["gene_pairs_ind_new"] = adata.uns["gene_pairs_ind_new"]
            old_adata.uns["gene_pairs_per_metabolite"] = adata.uns["gene_pairs_per_metabolite"]
            old_adata.uns["gene_pairs"] = adata.uns["gene_pairs"]
            old_adata.uns["gene_pairs_per_ct_pair"] = adata.uns["gene_pairs_per_ct_pair"]
            old_adata.uns["cell_type_pairs"] = adata.uns["cell_type_pairs"]
            old_adata.uns["cell_type_list"] = adata.uns["cell_type_list"]
            old_adata.uns["genes"] = adata.uns["genes"]
            if test == 'parametric':
                old_adata.uns['lc_zs_3d'] = adata.uns['lc_zs_3d']
            old_adata.uns['lcs_3d'] = adata.uns['lcs_3d']

            if return_all_adatas:
                return [old_adata, adata]
            else:
                return old_adata
        else:
            return adata

    # compute gene pair/metabolite scores

    # # For each metabolite and gene pair, an overall score is computed for every cell summarizing the presence of a given metabolite or the presence of interaction driven by two genes.
    # compute_communication_scores(
    #     adata,
    #     layer_key,
    #     deconv_data,
    #     cell_type_list,
    # )

    # # Computes local autocorrelation for every metabolite/gene pair. Assess the statistical significance.
    # compute_communication_autocorrelation(
    #     adata,
    #     compute_neighbors_on_key,
    # )

    # compute_communication_modules(
    #     adata,
    #     layer_key
    # )

    # data_accessor.compute_one_vs_all_scores() # Compute one vs all DE analysis of the metabolite scores for every categorical variable.
    # data_accessor.compute_gene_score_per_metabolite() # Compute gene score per metabolite. Computes the correlation between gene expression and metabolite scores for every gene-metabolite pair.
