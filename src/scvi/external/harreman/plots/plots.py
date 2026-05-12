import warnings
from typing import Optional, Union, Literal
from anndata import AnnData
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import leaves_list


def local_correlation_plot(
    adata: AnnData,
    mod_cmap='tab10',
    vmin=-10,
    vmax=10,
    z_cmap='RdBu_r',
    yticklabels=False,
    use_super_modules=False,
    show=True,
):
    """
    Plot a hierarchical-clustered heatmap of pairwise correlation Z-scores.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing:
        - ``uns['lc_zs']``: DataFrame of pairwise correlation Z-scores
        - ``uns['modules']`` or ``uns['super_modules']``
        - ``uns['linkage']``: scipy linkage matrix
    mod_cmap : str, default "tab10"
        Colormap for module color annotations.
    vmin, vmax : float, default -10 and 10, respectively
        Limits for the heatmap color scale.
    z_cmap : str, default "RdBu_r"
        Colormap for the Z-score heatmap.
    yticklabels : bool, default False
        Whether to show y-axis tick labels.
    use_super_modules : bool, default False
        Whether to use ``uns['super_modules']`` instead of ``uns['modules']``.
    show : bool, default True
        If ``True``, display the plot.
    """
    
    local_correlation_z = adata.uns["lc_zs"]
    modules = adata.uns["super_modules"] if use_super_modules else adata.uns["modules"]
    linkage = adata.uns["linkage"]

    row_colors = None
    colors = list(plt.get_cmap(mod_cmap).colors)
    module_colors = {i: colors[(i-1) % len(colors)] for i in modules.unique()}
    module_colors[-1] = '#ffffff'
    
    modules = modules[local_correlation_z.index]

    row_colors1 = pd.Series(
        [module_colors[i] for i in modules],
        index=local_correlation_z.index,
    )

    row_colors = pd.DataFrame({
        "Modules": row_colors1,
    })

    cm = sns.clustermap(
        local_correlation_z,
        row_linkage=linkage,
        col_linkage=linkage,
        vmin=vmin,
        vmax=vmax,
        cmap=z_cmap,
        xticklabels=False,
        yticklabels=yticklabels,
        row_colors=row_colors,
        rasterized=True,
    )

    fig = plt.gcf()
    plt.sca(cm.ax_heatmap)
    plt.ylabel("")
    plt.xlabel("")

    cm.ax_row_dendrogram.remove()

    # Add 'module X' annotations
    ii = leaves_list(linkage)

    mod_reordered = modules.iloc[ii]
    
    adata.uns['mod_reordered'] = [mod for mod in mod_reordered.unique() if mod != -1]

    mod_map = {}
    y = np.arange(modules.size)

    for x in mod_reordered.unique():
        if x == -1:
            continue

        mod_map[x] = y[mod_reordered == x].mean()

    plt.sca(cm.ax_row_colors)
    for mod, mod_y in mod_map.items():
        plt.text(-.5, y=mod_y, s="Module {}".format(mod),
                 horizontalalignment='right',
                 verticalalignment='center')
    plt.xticks([])

    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    min_aa.set_ylabel('Z-Scores')
    min_aa.yaxis.set_label_position("left")
    
    if show:
        plt.show()


def average_local_correlation_plot(
    adata: AnnData,
    mod_cmap='tab10',
    vmin=-10,
    vmax=10,
    cor_cmap='RdBu_r',
    yticklabels=False,
    row_cluster=True,
    col_cluster=True,
    use_super_modules=False,
    super_module_dict=None,
    show=True,
):
    """
    Plot the average pairwise correlation Z-scores between modules.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing:
        - ``uns['lc_zs']``: DataFrame of pairwise correlation Z-scores
        - ``uns['modules']`` or ``uns['super_modules']``  
        - ``uns['mod_reordered']``: ordering from `local_correlation_plot`
    mod_cmap : str, default "tab10"
        Colormap for module annotations.
    vmin, vmax : float, default -10 and 10, respectively
        Color scale limits for Z-scores.
    cor_cmap : str, default "RdBu_r"
        Colormap for the averaged correlation matrix.
    yticklabels : bool, default False
        Whether to display module labels on the heatmap.
    row_cluster, col_cluster : bool, default True
        Whether to apply clustering along rows/columns.
    use_super_modules : bool, default False
        Whether to use super-modules.
    super_module_dict : dict, optional
        Map of super-module → list of modules, used to color modules by their
        parent super-module.
    show : bool, default True
        Whether to display the plot.
    """
    
    local_correlation_z = adata.uns["lc_zs"]
    modules = adata.uns["super_modules"] if use_super_modules else adata.uns["modules"]
    
    avg_local_correlation_z = local_correlation_z.copy()
    avg_local_correlation_z['module_row'] = modules
    avg_local_correlation_z = avg_local_correlation_z.set_index('module_row', append=True)
    avg_local_correlation_z.columns = pd.MultiIndex.from_arrays([modules[avg_local_correlation_z.columns].values, avg_local_correlation_z.columns])

    avg_local_correlation_z = avg_local_correlation_z.groupby(level=1).mean().groupby(level=0, axis=1).mean()
    avg_local_correlation_z = avg_local_correlation_z.loc[avg_local_correlation_z.index != -1, avg_local_correlation_z.columns != -1]
    mod_reordered = adata.uns['mod_reordered']
    avg_local_correlation_z = avg_local_correlation_z.loc[mod_reordered, mod_reordered]

    row_colors = None
    colors = list(plt.get_cmap(mod_cmap).colors)
    if super_module_dict:
        module_colors = {
            mod: colors[(sm-1) % len(colors)]
            for sm, mods in super_module_dict.items()
            for mod in mods
        }
    else:
        module_colors = {mod: colors[(mod-1) % len(colors)] for mod in modules.unique()}
    module_colors[-1] = '#ffffff'

    row_colors = pd.DataFrame({
        "Modules": module_colors,
    })

    cm = sns.clustermap(
        avg_local_correlation_z,
        vmin=vmin,
        vmax=vmax,
        cmap=cor_cmap,
        xticklabels=False,
        yticklabels=yticklabels,
        row_colors=row_colors,
        rasterized=True,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
    )

    fig = plt.gcf()
    plt.sca(cm.ax_heatmap)
    plt.ylabel("")
    plt.xlabel("")

    cm.ax_row_dendrogram.remove()

    if row_cluster:
        reordered_indices = cm.dendrogram_row.reordered_ind
        mod_reordered = [avg_local_correlation_z.index[i] for i in reordered_indices]

    mod_map = {}
    y = np.arange(len(mod_reordered))

    for x in mod_reordered:
        if x == -1:
            continue

        mod_map[x] = y[mod_reordered == x].mean() + 0.5

    plt.sca(cm.ax_row_colors)
    for mod, mod_y in mod_map.items():
        plt.text(-.5, y=mod_y, s="Module {}".format(mod),
                    horizontalalignment='right',
                    verticalalignment='center')
    plt.xticks([])

    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    min_aa.set_ylabel('Avg. local correlation Z')
    min_aa.yaxis.set_label_position("left")

    if show:
        plt.show()


def module_score_correlation_plot(
    adata: AnnData,
    mod_cmap='tab10',
    vmin=-1,
    vmax=1,
    cor_cmap='RdBu_r',
    yticklabels=False,
    method='pearson',
    use_super_modules=False,
    super_module_dict=None,
    row_cluster=True,
    col_cluster=True,
    show=True,
):
    """
    Plot correlations between module scores across cells.

    Parameters
    ----------
    adata : AnnData
        Must contain:
        - ``obsm['module_scores']`` or ``obsm['super_module_scores']``
        - ``uns['modules']`` or ``uns['super_modules']``
        - ``uns['mod_reordered']``
    mod_cmap : str, default "tab10"
        Colormap for module annotations.
    vmin, vmax : float, default -1 and 1, respectively
        Color scale limits for the correlation heatmap.
    cor_cmap : str, default "RdBu_r"
        Colormap for correlation values.
    yticklabels : bool, default False
        Whether to show y-axis labels.
    method : {"pearson", "spearman"}, default "pearson"
        Correlation method.
    use_super_modules : bool, default False
        Whether to use super-module scores.
    super_module_dict : dict, optional
        Coloring scheme based on parent super-modules.
    row_cluster, col_cluster : bool, default True
        Whether to cluster rows/columns before plotting.
    show : bool, default True
        Whether to display the plot.
    """
    
    module_scores = adata.obsm['super_module_scores'] if use_super_modules else adata.obsm['module_scores']
    modules = adata.uns["super_modules"] if use_super_modules else adata.uns["modules"]

    cor_matrix = module_scores.corr(method)
    mod_int = [int(mod.split(' ')[1]) for mod in cor_matrix.index]
    cor_matrix.index = cor_matrix.columns = mod_int
    mod_reordered = adata.uns['mod_reordered']
    cor_matrix = cor_matrix.loc[mod_reordered, mod_reordered]

    row_colors = None
    colors = list(plt.get_cmap(mod_cmap).colors)
    if super_module_dict:
        module_colors = {
            mod: colors[(sm-1) % len(colors)]
            for sm, mods in super_module_dict.items()
            for mod in mods
        }
    else:
        module_colors = {mod: colors[(mod-1) % len(colors)] for mod in modules.unique()}
    module_colors[-1] = '#ffffff'

    row_colors = pd.DataFrame({
        "Modules": module_colors,
    })

    cm = sns.clustermap(
        cor_matrix,
        vmin=vmin,
        vmax=vmax,
        cmap=cor_cmap,
        xticklabels=False,
        yticklabels=yticklabels,
        row_colors=row_colors,
        rasterized=True,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
    )

    fig = plt.gcf()
    plt.sca(cm.ax_heatmap)
    plt.ylabel("")
    plt.xlabel("")

    cm.ax_row_dendrogram.remove()

    if row_cluster:
        reordered_indices = cm.dendrogram_row.reordered_ind
        mod_reordered = [cor_matrix.index[i] for i in reordered_indices]

    mod_map = {}
    y = np.arange(len(mod_reordered))

    for x in mod_reordered:
        if x == -1:
            continue

        mod_map[x] = y[mod_reordered == x].mean() + 0.5

    plt.sca(cm.ax_row_colors)
    for mod, mod_y in mod_map.items():
        plt.text(-.5, y=mod_y, s="Module {}".format(mod),
                    horizontalalignment='right',
                    verticalalignment='center')
    plt.xticks([])

    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    min_aa.set_ylabel(f'{method.capitalize()} R')
    min_aa.yaxis.set_label_position("left")
    
    if show:
        plt.show()


def plot_interacting_cell_scores(
    adata: AnnData,
    interactions: Optional[list] = None,
    coords_obsm_key: Optional[str] = None,
    test: Optional[Union[Literal["parametric"], Literal["non-parametric"]]] = None,
    only_sig_values: Optional[bool] = False,
    use_FDR: Optional[bool] = True,
    normalize_values: Optional[bool] = False,
    sample_specific: Optional[bool] = False,
    s: Optional[float] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Optional[tuple] = (10,10),
    cmap: Optional[str] = 'Reds',
    colorbar: Optional[bool] = True,
    swap_y_axis: Optional[bool] = False,
):
    """
    Plot spatial maps of interacting cell scores for selected gene pairs or metabolites.

    Parameters
    ----------
    adata : AnnData
        Must contain interacting cell statistics in:
        ``uns['interacting_cell_results'][test]['gp'/'m']``.
    interactions : list of str
        Gene pairs or metabolites to plot.
    coords_obsm_key : str
        Key in ``adata.obsm`` containing spatial coordinates.
    test : {"parametric", "non-parametric"}
        Determines which statistical results to load.
    only_sig_values : bool, default False
        If True, plot only significant values (FDR or p-value).
    use_FDR : bool, default True
        If ``only_sig_values=True``, choose FDR instead of raw p-values.
    normalize_values : bool, default False
        Apply per-interaction min–max normalization.
    sample_specific : bool, default False
        Plot each sample separately, using ``uns['sample_key']``.
    s : float
        Dot size.
    vmin, vmax : float or str ("p5", "p95"), optional
        Color scale limits; percentiles allowed.
    figsize: tuple, default (10,10)
        Figure size.
    cmap : str, default "Reds"
        Colormap for the score intensity.
    colorbar : bool, default True
        Show or hide the colorbar.
    swap_y_axis : bool, default False
        Flip the y-axis for visualization conventions.
    """
    
    if isinstance(vmin, str) and 'p' not in vmin:
        raise ValueError('"vmin" needs to be either a numeric value or a percentile: e.g. "p5".')
    if isinstance(vmax, str) and 'p' not in vmax:
        raise ValueError('"vmax" needs to be either a numeric value or a percentile: e.g. "p95".')
    
    if test not in ['parametric', 'non-parametric']:
        raise ValueError('The "test" variable should be one of ["parametric", "non-parametric"].')
    
    test_str = 'p' if test == 'parametric' else 'np'
    
    if sample_specific and 'sample_key' not in adata.uns.keys():
        raise ValueError('Sample information not found. Run Harreman using the "sample_key" parameter.')
    
    if only_sig_values:
        sig_str = 'FDR' if use_FDR else 'pval'
        interacting_cell_scores_gp = adata.uns['interacting_cell_results'][test_str]['gp'][f'cs_sig_{sig_str}']
        interacting_cell_scores_m = adata.uns['interacting_cell_results'][test_str]['m'][f'cs_sig_{sig_str}']
    else:
        interacting_cell_scores_gp = adata.uns['interacting_cell_results'][test_str]['gp']['cs']
        interacting_cell_scores_m = adata.uns['interacting_cell_results'][test_str]['m']['cs']
    
    if interactions is None:
        raise ValueError("Please provide a LR pair or a metabolite.")
    
    interacting_cell_scores_gp = pd.DataFrame(interacting_cell_scores_gp, index=adata.obs_names, columns=adata.uns['gene_pairs_sig_names'])
    interacting_cell_scores_m = pd.DataFrame(interacting_cell_scores_m, index=adata.obs_names, columns=adata.uns['metabolites'])
    gene_pairs = [inter for inter in interactions if inter in adata.uns['gene_pairs_sig_names']]
    metabs = [inter for inter in interactions if inter in adata.uns['metabolites']]
    if len(gene_pairs) > 0 and len(metabs) > 0:
        interacting_cell_scores = pd.concat([interacting_cell_scores_gp, interacting_cell_scores_m], axis=1)
    elif len(gene_pairs) == 0 and len(metabs) == 0:
        raise ValueError("The provided LR pairs and/or metabolites don't have significant interactions.")
    else:
        interacting_cell_scores = interacting_cell_scores_gp if len(gene_pairs) > 0 else interacting_cell_scores_m
    
    scores = interacting_cell_scores[interactions]
    
    if normalize_values:
        scores = scores.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0) #We apply min-max normalization
    
    for interaction in interactions:
        if isinstance(vmin, str):
            vmin_new = int(vmin.split('p')[1])
            vmin_new = np.percentile(scores[interaction], vmin_new)
        else:
            vmin_new = vmin
        if isinstance(vmax, str):
            vmax_new = int(vmax.split('p')[1])
            vmax_new = np.percentile(scores[interaction], vmax_new)
        else:
            vmax_new = vmax
        
        if sample_specific:
            sample_key = adata.uns['sample_key']
            for sample in adata.obs[sample_key].unique().tolist():
                print(sample)
                plot_interaction(adata[adata.obs[sample_key] == sample], scores.loc[adata.obs[sample_key] == sample], interaction, coords_obsm_key, s, vmin_new, vmax_new, figsize, cmap, colorbar, swap_y_axis)
                plt.show()
                plt.close()
        else:
            plot_interaction(adata, scores, interaction, coords_obsm_key, s, vmin_new, vmax_new, figsize, cmap, colorbar, swap_y_axis)
            plt.show()
            plt.close()


def plot_ct_interacting_cell_scores(
    adata: AnnData,
    deconv_adata: Optional[AnnData] = None,
    cell_type_pair: Optional[list] = None,
    interactions: Optional[list] = None,
    coords_obsm_key: Optional[str] = None,
    test: Optional[Union[Literal["parametric"], Literal["non-parametric"]]] = None,
    agg_only: Optional[bool] = False,
    normalize_values: Optional[bool] = False,
    sample_specific: Optional[bool] = False,
    s: Optional[float] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Optional[tuple] = (10,10),
    cmap: Optional[str] = 'Reds',
    colorbar: Optional[bool] = True,
    swap_y_axis: Optional[bool] = False,
):
    """
    Plot cell-type–specific interacting cell scores (for gene pairs or metabolites) across spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        AnnData containing per-cell-type interaction scores.
    deconv_adata : AnnData, optional
        If provided, interaction results are copied from this object into `adata`.
    cell_type_pair : list of str or list of tuple
        Cell type pairs to visualize  
        - tuple: ("T cell", "Macrophage")  
        - string: "T cell" (matches any pair containing this cell type)
    interactions : list of str
        Gene pairs or metabolites to visualize.
    coords_obsm_key : str
        Key for spatial coordinates in ``obsm``.
    test : {"parametric", "non-parametric"}
        Statistical test used during computation.
    agg_only : bool, default False
        Whether to plot only aggregated per-cell-type interactions.
    normalize_values : bool, default False
        Apply min–max normalization to each column.
    sample_specific : bool, default False
        Plot one figure per sample.
    s : float
        Dot size.
    vmin, vmax : float or str ("p5"), optional
        Color scale limits; percentiles allowed.
    figsize : tuple, default (10,10)
        Figure size.
    cmap : str, default "Reds"
        Colormap for plotting.
    colorbar : bool, default True
        Whether to display the colorbar.
    swap_y_axis : bool, default False
        Flip y-axis orientation.
    """
    
    if isinstance(vmin, str) and 'p' not in vmin:
        raise ValueError('"vmin" needs to be either a numeric value or a percentile: e.g. "p5".')
    if isinstance(vmax, str) and 'p' not in vmax:
        raise ValueError('"vmax" needs to be either a numeric value or a percentile: e.g. "p95".')
    
    if test not in ['parametric', 'non-parametric']:
        raise ValueError('The "test" variable should be one of ["parametric", "non-parametric"].')
    
    test_str = 'p' if test == 'parametric' else 'np'
    
    if sample_specific and 'sample_key' not in adata.uns.keys():
        raise ValueError('Sample information not found. Run Harreman using the "sample_key" parameter.')
    
    if deconv_adata is not None:
        adata.uns[f'ct_interacting_cell_results_{test_str}_gp_cs_df'] = deconv_adata.uns[f'ct_interacting_cell_results_{test_str}_gp_cs_df']
        adata.uns[f'ct_interacting_cell_results_{test_str}_m_cs_df'] = deconv_adata.uns[f'ct_interacting_cell_results_{test_str}_m_cs_df']
    
    if interactions is None:
        raise ValueError("Please provide a LR pair or a metabolite.")
    
    cell_type_pair = [] if cell_type_pair is None else cell_type_pair

    if not isinstance(cell_type_pair, list):
        raise ValueError(
            'The "cell_type_pair" variable must be None, a list of strings, or a list of tuples.'
        )
        
    ct_pairs = []
    if cell_type_pair:
        for ct in cell_type_pair:
            if isinstance(ct, tuple):
                ct_pairs.append(f"{ct[0]} - {ct[1]}")
            elif isinstance(ct, str):
                ct_pairs.append(ct)
            else:
                raise ValueError(
                    'Each element in "cell_type_pair" must be either a tuple or a string.'
                )
    
    interacting_cell_scores_gp = adata.obsm[f'ct_interacting_cell_results_{test_str}_gp_cs_df']
    interacting_cell_scores_m = adata.obsm[f'ct_interacting_cell_results_{test_str}_m_cs_df']
    
    def match_columns(df, ct_pairs, interactions):
        matched_columns = []
        per_cell_aggregation = {}
        for col in df.columns:
            ct_pair_str, interaction = col.split(': ', 1)
            if interaction not in interactions:
                continue
            if not ct_pairs:
                matched_columns.append(col)
                continue
            try:
                ct1, ct2 = ct_pair_str.split(' - ')
            except ValueError:
                continue  # skip malformed cell type pairs

            for query in ct_pairs:
                if ' - ' in query:
                    # Exact match of cell type pair
                    if query == f"{ct1} - {ct2}":
                        matched_columns.append(col)
                        break
                else:
                    # Single cell type: match if in either position
                    if query == ct1 or query == ct2:
                        matched_columns.append(col)
                        key = f"{query}: {interaction}"
                        if key not in per_cell_aggregation:
                            per_cell_aggregation[key] = df[col].copy()
                        else:
                            per_cell_aggregation[key] += df[col]
                        break

        selected_df = df[matched_columns]
        if per_cell_aggregation:
            agg_df = pd.concat(per_cell_aggregation, axis=1)
        else:
            agg_df = pd.DataFrame(index=df.index)  # empty fallback

        return selected_df, agg_df
    
    gp_selected, gp_aggregated = match_columns(interacting_cell_scores_gp, ct_pairs, interactions)
    m_selected, m_aggregated = match_columns(interacting_cell_scores_m, ct_pairs, interactions)
    
    dfs = [gp_aggregated, m_aggregated] if agg_only else [gp_selected, gp_aggregated, m_selected, m_aggregated]
    scores = pd.concat(dfs, axis=1)
    
    if normalize_values:
        scores = scores.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0) #We apply min-max normalization
    
    for interaction in scores.columns:
        if isinstance(vmin, str):
            vmin_new = float(vmin.split('p')[1])
            vmin_new = np.percentile(scores[interaction], vmin_new)
        else:
            vmin_new = vmin
        if isinstance(vmax, str):
            vmax_new = float(vmax.split('p')[1])
            vmax_new = np.percentile(scores[interaction], vmax_new)
        else:
            vmax_new = vmax
        
        if sample_specific:
            sample_key = adata.uns['sample_key']
            for sample in adata.obs[sample_key].unique().tolist():
                print(sample)
                plot_ct_interaction(adata[adata.obs[sample_key] == sample], scores.loc[adata.obs[sample_key] == sample], interaction, coords_obsm_key, s, vmin_new, vmax_new, figsize, cmap, colorbar, swap_y_axis)
                plt.show()
                plt.close()
        else:
            plot_ct_interaction(adata, scores, interaction, coords_obsm_key, s, vmin_new, vmax_new, figsize, cmap, colorbar, swap_y_axis)
            plt.show()
            plt.close()


def plot_interaction(adata, scores, interaction, coords_obsm_key, s, vmin, vmax, figsize, cmap, colorbar, swap_y_axis):

    if isinstance(adata.obsm[coords_obsm_key], pd.DataFrame):
        coords = adata.obsm[coords_obsm_key].values
    else:
        coords = adata.obsm[coords_obsm_key]
    
    ax = plt.subplot(111)
    ax.set_aspect('equal', adjustable='box')
    _prettify_axis(ax, spatial=True)
    if swap_y_axis:
        plt.scatter(coords[:,0], -coords[:,1], c=scores[interaction], cmap=cmap, s=s, vmin=vmin, vmax=vmax)
    else:
        plt.scatter(coords[:,0], coords[:,1], c=scores[interaction], cmap=cmap, s=s, vmin=vmin, vmax=vmax)
    plt.title(interaction)
    if colorbar:
        plt.colorbar()
        

def plot_ct_interaction(adata, scores, interaction, coords_obsm_key, s, vmin, vmax, figsize, cmap, colorbar, swap_y_axis):

    if isinstance(adata.obsm[coords_obsm_key], pd.DataFrame):
        coords = adata.obsm[coords_obsm_key].values
    else:
        coords = adata.obsm[coords_obsm_key]
    
    ax = plt.subplot(111)
    ax.set_aspect('equal', adjustable='box')
    _prettify_axis(ax, spatial=True)
    if swap_y_axis:
        plt.scatter(coords[:,0], -coords[:,1], c=scores[interaction], cmap=cmap, s=s, vmin=vmin, vmax=vmax)
    else:
        plt.scatter(coords[:,0], coords[:,1], c=scores[interaction], cmap=cmap, s=s, vmin=vmin, vmax=vmax)
    plt.title(interaction)
    if colorbar:
        plt.colorbar()


def _prettify_axis(ax, spatial=False):
    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    if spatial:
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("Spatial1")
        plt.ylabel("Spatial2")


def plot_signature_for_selection(adata, signature, coords_obsm_key, s, vmin, vmax, figsize, cmap, colorbar):

    scores = adata.obsm['vision_signatures']

    if isinstance(adata.obsm[coords_obsm_key], pd.DataFrame):
        coords = adata.obsm[coords_obsm_key].values
    else:
        coords = adata.obsm[coords_obsm_key]
    
    points = np.column_stack([coords[:,0], coords[:,1]])
    
    plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    _prettify_axis(ax, spatial=True)
    p = plt.scatter(coords[:,0], coords[:,1], c=scores[signature], cmap=cmap, s=s, vmin=vmin, vmax=vmax)
    plt.title(signature)
    if colorbar:
        plt.colorbar()
    
    return p, ax, points


def plot_selection_histplot(adata, signature, group):

    adata.obs['selected'] = group
    adata.obs['selected'][adata.obs['selected'] == 1] = 'Selection'
    adata.obs['selected'][adata.obs['selected'] == 0] = 'Remainder'

    if signature not in adata.obs:
        adata.obs[signature] = adata.obsm['vision_signatures'][signature]

    ax = sns.histplot(data=adata.obs, x=signature, hue=adata.obs['selected'].tolist(), bins=30, palette={'Selection': '#FF7F00', 'Remainder': '#1F78B4'})
    plt.show()

    return


def plot_vision_autocorrelation(
    adata, 
    type: Optional[Union[Literal["observations"], Literal["signatures"]]] = None,
    center: Optional[int] = 0.5,
    figsize: Optional[tuple] = (1,10),
    cmap: Optional[str] = 'coolwarm',
    cbar: Optional[bool] = True
):

    if type not in ["observations", "signatures"]:
        raise ValueError('The "type" variable should be one of ["observations", "signatures"].')
    
    type_str = 'vision_obs_df_scores' if type == "observations" else 'vision_signature_scores'

    masked_data = adata.uns[type_str][['c_prime']].where((adata.uns[type_str][['fdr']] < 0.05).values)
    masked_data = masked_data.sort_values('c_prime', ascending=False)
    masked_data.columns = ['Consistency']

    plt.figure(figsize=figsize)
    sns.heatmap(masked_data, annot=masked_data, cmap=cmap, fmt='.2f', cbar=cbar, center=center)
    plt.show()

    return


def plot_vision_de_results(
    adata,
    type: Optional[Union[Literal["observations"], Literal["signatures"]]] = None,
    var: str = None,
    center: Optional[int] = 0.5,
    figsize: Optional[tuple] = (3,10),
    cmap: Optional[str] = 'coolwarm',
    cbar: Optional[bool] = True
):

    if var is None:
        raise ValueError('The "var" variable should be a categorical variable to plot.')
    
    if type not in ["observations", "signatures"]:
        raise ValueError('The "type" variable should be one of ["observations", "signatures"].')
    
    type_score_str = f'one_vs_all_obs_cols_{var}_scores' if type == "observations" else f'one_vs_all_signatures_{var}_scores'
    type_pval_str = f'one_vs_all_obs_cols_{var}_pvals' if type == "observations" else f'one_vs_all_signatures_{var}_padj'

    mask = adata.uns[type_pval_str] < 0.05

    plt.figure(figsize=figsize)
    sns.heatmap(adata.uns[type_score_str], mask=~mask, cmap=cmap, annot=mask.applymap(lambda x: '*' if x else ''), fmt='', cbar=cbar, center=center)
    plt.show()

    return


def plot_sig_mod_correlation(
    adata,
    x_rotation: Optional[int] = 0,
    y_rotation: Optional[int] = 0,
    use_FDR: Optional[bool] = True,
    subset_signatures: Optional[list] = None,
    subset_modules: Optional[list] = None,
    cmap: Optional[str] = 'RdBu_r',
):
    
    coef = adata.uns['sig_mod_correlation_coefs'] if 'sig_mod_correlation_coefs' in adata.uns.keys() else None
    if use_FDR:
        padj = adata.uns['sig_mod_correlation_FDR'] if 'sig_mod_correlation_FDR' in adata.uns.keys() else None
    else:
        padj = adata.uns['sig_mod_correlation_pvals'] if 'sig_mod_correlation_pvals' in adata.uns.keys() else None

    if coef is None or padj is None:
        raise ValueError('Run the "harreman.hs.integrate_vision_hotspot_results" function before plotting the results.')
    
    coef = coef.loc[subset_signatures] if subset_signatures is not None else coef
    padj = padj.loc[subset_signatures] if subset_signatures is not None else padj
    coef = coef[subset_modules] if subset_modules is not None else coef
    padj = padj[subset_modules] if subset_modules is not None else padj
    
    coef = coef[padj < 0.05].dropna(how='all').copy()
    padj = padj[padj < 0.05].dropna(how='all').copy()

    coef.replace(np.nan, 0, inplace=True)
    padj.replace(np.nan, 1, inplace=True)
    
    cmap = mpl.colormaps.get_cmap(cmap)
    cmap.set_bad("gray")

    g = sns.clustermap(coef, cmap=cmap, xticklabels=True, yticklabels=True, mask=padj > 0.05, center=0)

    fig = plt.gcf()

    for tick in g.ax_heatmap.get_xticklabels():
        tick.set_rotation(x_rotation)
    for tick in g.ax_heatmap.get_yticklabels():
        tick.set_rotation(y_rotation)

    padj = padj[g.data2d.columns]

    for i, ix in enumerate(g.dendrogram_row.reordered_ind):
        for j in range(len(coef.columns)):
                text = g.ax_heatmap.text(
                    j + 0.5,
                    i + 0.5,
                    "***" if padj.iloc[ix, j] < 0.0005 else "**"
                    if padj.iloc[ix, j] < 0.005 else "*" if padj.iloc[ix, j] < 0.05 else '',
                    ha="center",
                    va="center",
                    color="black",
                )
    
    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    label = 'Spearman R' if adata.uns['cor_method'] == 'spearman' else 'Pearson R'
    min_aa.set_ylabel(label)
    min_aa.yaxis.set_label_position("left")
    
    plt.show()


def plot_sig_mod_enrichment(
    adata,
    x_rotation: Optional[int] = 0,
    y_rotation: Optional[int] = 0,
    use_FDR: Optional[bool] = True,
    subset_signatures: Optional[list] = None,
    subset_modules: Optional[list] = None,
    cmap: Optional[str] = 'RdBu_r',
):
    
    coef = adata.uns['sig_mod_enrichment_stats'] if 'sig_mod_enrichment_stats' in adata.uns.keys() else None
    if use_FDR:
        padj = adata.uns['sig_mod_enrichment_FDR'] if 'sig_mod_enrichment_FDR' in adata.uns.keys() else None
    else:
        padj = adata.uns['sig_mod_enrichment_pvals'] if 'sig_mod_enrichment_pvals' in adata.uns.keys() else None

    if coef is None or padj is None:
        raise ValueError('Run the "harreman.hs.integrate_vision_hotspot_results" function before plotting the results.')
    
    coef = coef.loc[subset_signatures] if subset_signatures is not None else coef
    padj = padj.loc[subset_signatures] if subset_signatures is not None else padj
    coef = coef[subset_modules] if subset_modules is not None else coef
    padj = padj[subset_modules] if subset_modules is not None else padj
    
    coef = coef[padj < 0.05].T.dropna(how='all').copy()
    padj = padj[padj < 0.05].T.dropna(how='all').copy()

    coef.replace(np.nan, 0, inplace=True)
    padj.replace(np.nan, 1, inplace=True)
    
    cmap = mpl.colormaps.get_cmap(cmap)
    cmap.set_bad("gray")

    g = sns.clustermap(coef, cmap=cmap, xticklabels=True, yticklabels=True, mask=padj > 0.05, center=0)

    for tick in g.ax_heatmap.get_xticklabels():
        tick.set_rotation(x_rotation)
    for tick in g.ax_heatmap.get_yticklabels():
        tick.set_rotation(y_rotation)

    padj = padj[g.data2d.columns]

    for i, ix in enumerate(g.dendrogram_row.reordered_ind):
        for j in range(len(coef.columns)):
                text = g.ax_heatmap.text(
                    j + 0.5,
                    i + 0.5,
                    "***" if padj.iloc[ix, j] < 0.0005 else "**"
                    if padj.iloc[ix, j] < 0.005 else "*" if padj.iloc[ix, j] < 0.05 else '',
                    ha="center",
                    va="center",
                    color="black",
                )


def plot_interaction_module_correlation(
    adata,
    x_rotation: Optional[int] = 0,
    y_rotation: Optional[int] = 0,
    use_FDR: Optional[bool] = True,
    subset_interactions: Optional[list] = None,
    subset_modules: Optional[list] = None,
    cmap: Optional[str] = 'RdBu_r',
    figsize: Optional[tuple] = (10,10),
    threshold: Optional[float] = None,
):
    
    coef = adata.uns['interaction_module_correlation_coefs'].T if 'interaction_module_correlation_coefs' in adata.uns.keys() else None
    if use_FDR:
        padj = adata.uns['interaction_module_correlation_FDR'].T if 'interaction_module_correlation_FDR' in adata.uns.keys() else None
    else:
        padj = adata.uns['interaction_module_correlation_pvals'].T if 'interaction_module_correlation_pvals' in adata.uns.keys() else None

    if coef is None or padj is None:
        raise ValueError('Run the "harreman.tl.compute_interaction_module_correlation" function before plotting the results.')
    
    coef = coef.loc[subset_interactions] if subset_interactions is not None else coef
    padj = padj.loc[subset_interactions] if subset_interactions is not None else padj
    coef = coef[subset_modules] if subset_modules is not None else coef
    padj = padj[subset_modules] if subset_modules is not None else padj
    
    coef = coef[padj < 0.05].dropna(how='all').copy()
    padj = padj[padj < 0.05].dropna(how='all').copy()

    coef.replace(np.nan, 0, inplace=True)
    padj.replace(np.nan, 1, inplace=True)
    
    if threshold:
        padj = padj[(coef > threshold).any(axis=1)]
        coef = coef[(coef > threshold).any(axis=1)]
    
    cmap = mpl.colormaps.get_cmap(cmap)
    cmap.set_bad("gray")

    g = sns.clustermap(coef, cmap=cmap, xticklabels=True, yticklabels=True, mask=padj > 0.05, center=0, figsize=figsize)

    fig = plt.gcf()

    for tick in g.ax_heatmap.get_xticklabels():
        tick.set_rotation(x_rotation)
    for tick in g.ax_heatmap.get_yticklabels():
        tick.set_rotation(y_rotation)

    padj = padj[g.data2d.columns]

    for i, ix in enumerate(g.dendrogram_row.reordered_ind):
        for j in range(len(coef.columns)):
            text = g.ax_heatmap.text(
                j + 0.5,
                i + 0.5,
                "***" if padj.iloc[ix, j] < 0.0005 else "**"
                if padj.iloc[ix, j] < 0.005 else "*" if padj.iloc[ix, j] < 0.05 else '',
                ha="center",
                va="center",
                color="black",
            )
    
    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    label = 'Spearman R' if adata.uns['cor_method'] == 'spearman' else 'Pearson R'
    min_aa.set_ylabel(label)
    min_aa.yaxis.set_label_position("left")
    
    plt.show()
