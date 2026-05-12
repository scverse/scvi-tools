from typing import Optional
import time
import anndata
import numpy as np
import pandas as pd
from scipy.sparse import issparse


def counts_from_anndata(adata, layer_key=None, dense=False):
    # 1. Extract matrix
    if layer_key is None:
        counts = adata.X
    elif layer_key == "use_raw":
        counts = adata.raw.X
    else:
        counts = adata.layers[layer_key]

    # 2. Transpose efficiently
    if issparse(counts):
        counts = counts.transpose().tocsr(copy=False)  # keep CSR format for efficient row slicing
    else:
        counts = counts.T  # transpose numpy array directly

    # 3. Convert to dense if requested
    if dense:
        if issparse(counts):
            counts = counts.toarray()
        else:
            counts = np.asarray(counts)

    return counts


def setup_anndata(
        input_adata: anndata.AnnData,
        cell_types: list,
        compute_neighbors_on_key: str,
        cell_type_key: str,
        database_varm_key: str,
        sample_key: Optional[str],
        spot_diameter: int,
) -> anndata.AnnData:
    
    barcode_key = 'barcodes'
    database = input_adata.varm[database_varm_key]
    obs_names = input_adata.obs_names
    var_names = input_adata.var_names
    obsm_neighbors = input_adata.obsm[compute_neighbors_on_key]
    
    adatas = []
    for ct in cell_types:
        if ct not in input_adata.layers:
            continue
        X = input_adata.layers[ct]
        obs = {
            barcode_key: obs_names.values,
            cell_type_key: [ct] * len(obs_names)
        }
        if sample_key is not None:
            obs[sample_key] = input_adata.obs[sample_key].values

        adata = anndata.AnnData(
            X=X,
            obs=obs,
            var=pd.DataFrame(index=var_names),
            obsm={compute_neighbors_on_key: obsm_neighbors},
        )
        adata.obs_names = [f"{name}_{ct}" for name in obs_names]
        adatas.append(adata)

    out_adata = anndata.concat(adatas, axis=0)

    out_adata.uns.update({
        'database_varm_key': database_varm_key,
        'spot_diameter': spot_diameter,
        'barcode_key': barcode_key
    })
    out_adata.varm[database_varm_key] = database

    # Remove empty rows efficiently
    nonzero_mask = out_adata.X.sum(axis=1).A1 > 0 if hasattr(out_adata.X, 'A1') else out_adata.X.sum(axis=1) > 0
    out_adata._inplace_subset_obs(nonzero_mask)

    return out_adata


def setup_deconv_adata(
    adata: anndata.AnnData,
    compute_neighbors_on_key: Optional[str] = None,
    sample_key: Optional[str] = None,
    cell_type_list: Optional[list] = None,
    cell_type_key: Optional[str] = None,
    spot_diameter: Optional[int] = None,
    verbose: Optional[bool] = False,
):
    """Set up deconvolution AnnData.

    Parameters
    ----------
    adata
        AnnData object.
    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use neighbors stored in `adata`. If no neighbors have been previously computed an error will be raised.
    sample_key
        Sample information in case the data contains different samples or samples from different conditions. Input is key in `adata.obs`.
    cell_type_list
        Cell type or cluster information for the cell-cell communication analysis. Input is a list of keys in `adata.layers`.
    cell_type_key
        Cell type or cluster information for the cell-cell communication analysis. Input is key in `adata.obs`.
    spot_diameter
        Spot diameter of the spatial technology the dataset belongs to.
    verbose : bool, optional (default: False)
        Whether to print progress and status messages.

    """

    start = time.time()
    if verbose:
        print("Setting up deconvolution data...")

    uns = adata.uns
    deconv_adata = setup_anndata(
        input_adata=adata,
        cell_types=cell_type_list,
        compute_neighbors_on_key=compute_neighbors_on_key,
        cell_type_key=cell_type_key,
        database_varm_key=uns['database_varm_key'],
        sample_key=sample_key,
        spot_diameter=spot_diameter,
    )
    
    deconv_adata.uns.update({
        'cell_type_key': cell_type_key,
        'layer_key': None,
        'deconv_data': True,
        'database': uns["database"]
    })
        
    if uns["database"] in {'LR', 'both'}:
        deconv_adata.uns.update({
            'ligand': uns['ligand'],
            'receptor': uns['receptor'],
            'LR_database': uns['LR_database']
        })

    if uns["database"] in {'transporter', 'both'}:
        deconv_adata.uns.update({
            'importer': uns['importer'],
            'exporter': uns['exporter'],
            'import_export': uns['import_export'],
            'num_metabolites': uns['num_metabolites'],
            'metabolite_database': uns['metabolite_database']
        })
    
    if verbose:
        print("Finished setting up deconvolution data in %.3f seconds" %(time.time()-start))

    return adata, deconv_adata


def modify_uns_hotspot(adata):
    if 'modules' in adata.uns.keys():
        adata.var['modules'] = adata.uns['modules']
        del adata.uns['modules']
    
    if 'super_modules' in adata.uns.keys():
        adata.var['super_modules'] = adata.uns['super_modules']
        del adata.uns['super_modules']
    
    if 'lc_zs' in adata.uns.keys():
        genes = [' - '.join(gene) if isinstance(gene, tuple) else gene for gene in adata.uns['lc_zs'].columns]
        adata.uns['lc_zs'].index = genes
        adata.uns['lc_zs'].columns = genes
    
    return


def modify_uns_harreman(adata):
    uns_keys = ['ligand', 'receptor', 'LR_database', 'import_export']
    for uns_key in uns_keys:
        if uns_key in adata.uns.keys():
            adata.uns[uns_key] = adata.uns[uns_key].fillna('NA')
    
    if 'LR_database' in adata.uns.keys():
        adata.uns['LR_database'].columns = [col.replace('.', '_') for col in adata.uns['LR_database'].columns]
        adata.uns['LR_database']['ligand_transmembrane'] = adata.uns['LR_database']['ligand_transmembrane'].astype(str)
        adata.uns['LR_database']['receptor_transmembrane'] = adata.uns['LR_database']['receptor_transmembrane'].astype(str)

    if 'gene_pairs' in adata.uns.keys():
        gene_pairs_tmp = [(x, ' - '.join(y) if isinstance(y, (list, tuple)) else y) for x, y in adata.uns['gene_pairs']]
        gene_pairs_tmp = [(' - '.join(x) if isinstance(x, (list, tuple)) else x, y) for x, y in gene_pairs_tmp]
        adata.uns['gene_pairs'] = ['_'.join(gp) for gp in gene_pairs_tmp]
    
    if 'gene_pairs_sig' in adata.uns.keys():
        gene_pairs_sig_tmp = [(x, ' - '.join(y) if isinstance(y, (list, tuple)) else y) for x, y in adata.uns['gene_pairs_sig']]
        gene_pairs_sig_tmp = [(' - '.join(x) if isinstance(x, (list, tuple)) else x, y) for x, y in gene_pairs_sig_tmp]
        adata.uns['gene_pairs_sig'] = ['_'.join(gp) for gp in gene_pairs_sig_tmp]

    if 'gene_pairs_ind' in adata.uns.keys():
        gene_pairs_ind_tmp = [(x, ' - '.join(map(str, y)) if isinstance(y, (list, tuple)) else str(y)) for x, y in adata.uns['gene_pairs_ind']]
        gene_pairs_ind_tmp = [(' - '.join(map(str, x)) if isinstance(x, (list, tuple)) else str(x), y) for x, y in gene_pairs_ind_tmp]
        adata.uns['gene_pairs_ind'] = ['_'.join(gp) for gp in gene_pairs_ind_tmp]
    
    if 'gene_pairs_sig_ind' in adata.uns.keys():
        gene_pairs_sig_ind_tmp = [(x, ' - '.join(map(str, y)) if isinstance(y, (list, tuple)) else str(y)) for x, y in adata.uns['gene_pairs_sig_ind']]
        gene_pairs_sig_ind_tmp = [(' - '.join(map(str, x)) if isinstance(x, (list, tuple)) else str(x), y) for x, y in gene_pairs_sig_ind_tmp]
        adata.uns['gene_pairs_sig_ind'] = ['_'.join(gp) for gp in gene_pairs_sig_ind_tmp]
    
    if 'gene_pairs_per_metabolite' in adata.uns.keys():
        adata.uns['gene_pairs_per_metabolite'] = {key: {
            'gene_pair': [(' - '.join(gp_1) if isinstance(gp_1, (list, tuple)) else gp_1, gp_2) for gp_1, gp_2 in subdict['gene_pair']],
            'gene_type': subdict['gene_type']
        } for key, subdict in adata.uns['gene_pairs_per_metabolite'].items()}
        adata.uns['gene_pairs_per_metabolite'] = {key: {
            'gene_pair': [(gp_1, ' - '.join(gp_2) if isinstance(gp_2, (list, tuple)) else gp_2) for gp_1, gp_2 in subdict['gene_pair']],
            'gene_type': subdict['gene_type']
        } for key, subdict in adata.uns['gene_pairs_per_metabolite'].items()}
    
    if 'gene_pairs_per_ct_pair' in adata.uns.keys():
        adata.uns['gene_pairs_per_ct_pair'] = {key: [(x, ' - '.join(y) if isinstance(y, (list, tuple)) else y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_per_ct_pair'].items()}
        adata.uns['gene_pairs_per_ct_pair'] = {key: [(' - '.join(x) if isinstance(x, (list, tuple)) else x, y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_per_ct_pair'].items()}
        adata.uns['gene_pairs_per_ct_pair'] = {' - '.join(key): value for key, value in adata.uns['gene_pairs_per_ct_pair'].items()}
    
    if 'gene_pairs_per_ct_pair_ind' in adata.uns.keys():
        adata.uns['gene_pairs_per_ct_pair_ind'] = {' - '.join(key): value for key, value in adata.uns['gene_pairs_per_ct_pair_ind'].items()}
    
    if 'gene_pairs_per_ct_pair_sig_ind' in adata.uns.keys():
        adata.uns['gene_pairs_per_ct_pair_sig_ind'] = {' - '.join(key): value for key, value in adata.uns['gene_pairs_per_ct_pair_sig_ind'].items()}
    
    if 'gene_pairs_ind_per_ct_pair' in adata.uns.keys():
        adata.uns['gene_pairs_ind_per_ct_pair'] = {key: [(x, ' - '.join(map(str, y)) if isinstance(y, (list, tuple)) else y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_ind_per_ct_pair'].items()}
        adata.uns['gene_pairs_ind_per_ct_pair'] = {key: [(' - '.join(map(str, x)) if isinstance(x, (list, tuple)) else x, y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_ind_per_ct_pair'].items()}
        adata.uns['gene_pairs_ind_per_ct_pair'] = {' - '.join(key): value for key, value in adata.uns['gene_pairs_ind_per_ct_pair'].items()}
    
    if 'gene_pairs_ind_per_ct_pair_sig' in adata.uns.keys():
        adata.uns['gene_pairs_ind_per_ct_pair_sig'] = {key: [(x, ' - '.join(map(str, y)) if isinstance(y, (list, tuple)) else y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_ind_per_ct_pair_sig'].items()}
        adata.uns['gene_pairs_ind_per_ct_pair_sig'] = {key: [(' - '.join(map(str, x)) if isinstance(x, (list, tuple)) else x, y) for x, y in tuples_list] for key, tuples_list in adata.uns['gene_pairs_ind_per_ct_pair_sig'].items()}
        adata.uns['gene_pairs_ind_per_ct_pair_sig'] = {' - '.join(key): value for key, value in adata.uns['gene_pairs_ind_per_ct_pair_sig'].items()}

    if 'ccc_results' in adata.uns.keys():
        adata.uns['ccc_results']['cell_com_df_gp'] = adata.uns['ccc_results']['cell_com_df_gp'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
        adata.uns['ccc_results']['cell_com_df_m'] = adata.uns['ccc_results']['cell_com_df_m'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
        if 'cell_com_df_gp_sig' in adata.uns['ccc_results'].keys():
            adata.uns['ccc_results']['cell_com_df_gp_sig'] = adata.uns['ccc_results']['cell_com_df_gp_sig'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
            adata.uns['ccc_results']['cell_com_df_m_sig'] = adata.uns['ccc_results']['cell_com_df_m_sig'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
    
    if 'ct_ccc_results' in adata.uns.keys():
        adata.uns['ct_ccc_results']['cell_com_df_gp'] = adata.uns['ct_ccc_results']['cell_com_df_gp'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
        adata.uns['ct_ccc_results']['cell_com_df_m'] = adata.uns['ct_ccc_results']['cell_com_df_m'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
        if 'cell_com_df_gp_sig' in adata.uns['ct_ccc_results'].keys():
            adata.uns['ct_ccc_results']['cell_com_df_gp_sig'] = adata.uns['ct_ccc_results']['cell_com_df_gp_sig'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)
            adata.uns['ct_ccc_results']['cell_com_df_m_sig'] = adata.uns['ct_ccc_results']['cell_com_df_m_sig'].applymap(lambda x: ' - '.join(x) if isinstance(x, (list, tuple)) else x)

    return


def write_h5ad(
    adata: anndata.AnnData, 
    filename: Optional[str] = None,
):
    """
    Save an AnnData object to disk with Harreman-compatible preprocessing.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to save.
    filename : str, optional
        Path to the output `.h5ad` file. If omitted, an error is raised.

    Notes
    -----
    The wrapper ensures that custom Harreman/Hotspot fields remain fully restorable
    when loading the file with `read_h5ad()`.
    """
    
    if filename is None:
        raise ValueError('Please provide the path to save the file.')
    elif not filename.endswith('h5ad'):
        filename = filename+'.h5ad'
    
    if 'distances' in adata.obsp.keys():
        adata.obsp['distances'] = adata.obsp['distances'].tocsr()
    
    modify_uns_hotspot(adata)
    modify_uns_harreman(adata)
    adata.write(filename)


def recover_uns_hotspot(adata):
    if 'modules' not in adata.uns and 'modules' in adata.var.columns:
        adata.uns['modules'] = adata.var['modules'].dropna().astype(int).copy()
        del adata.var['modules']

    if 'super_modules' not in adata.uns and 'super_modules' in adata.var.columns:
        adata.uns['super_modules'] = adata.var['super_modules'].dropna().astype(int).copy()
        del adata.var['super_modules']

    if 'lc_zs' in adata.uns:
        adata.uns['lc_zs'].index = [tuple(g.split(' - ')) if ' - ' in g else g for g in adata.uns['lc_zs'].index]
        adata.uns['lc_zs'].columns = [tuple(g.split(' - ')) if ' - ' in g else g for g in adata.uns['lc_zs'].columns]


def recover_uns_harreman(adata):
    uns_keys = ['ligand', 'receptor', 'LR_database', 'import_export']
    for uns_key in uns_keys:
        if uns_key in adata.uns.keys():
            adata.uns[uns_key] = adata.uns[uns_key].replace("NA", np.nan)
            
    if 'LR_database' in adata.uns:
        original_columns = ['interaction_name', 'pathway_name', 'agonist',
                            'antagonist', 'co_A_receptor', 'co_I_receptor', 'evidence',
                            'annotation', 'interaction_name_2', 'is_neurotransmitter',
                            'ligand.symbol', 'ligand.family', 'ligand.location', 'ligand.keyword',
                            'ligand.secreted_type', 'ligand.transmembrane', 'receptor.symbol',
                            'receptor.family', 'receptor.location', 'receptor.keyword',
                            'receptor.surfaceome_main', 'receptor.surfaceome_sub',
                            'receptor.adhesome', 'receptor.secreted_type', 'receptor.transmembrane',
                            'version']
        mod_columns = [col.replace('.', '_') for col in original_columns]
        adata.uns['LR_database'][mod_columns].columns = original_columns
        
    if 'gene_pairs' in adata.uns:
        gene_pairs_tmp = [tuple(gp.split('_')) for gp in adata.uns['gene_pairs']]
        gene_pairs_tmp = [(x, list(y.split(' - ')) if ' - ' in y else y) for x, y in gene_pairs_tmp]
        adata.uns['gene_pairs'] = [(list(x.split(' - ')) if ' - ' in x else x, y) for x, y in gene_pairs_tmp]
    
    if 'gene_pairs_sig' in adata.uns:
        gene_pairs_sig_tmp = [tuple(gp.split('_')) for gp in adata.uns['gene_pairs_sig']]
        gene_pairs_sig_tmp = [(x, list(y.split(' - ')) if ' - ' in y else y) for x, y in gene_pairs_sig_tmp]
        adata.uns['gene_pairs_sig'] = [(list(x.split(' - ')) if ' - ' in x else x, y) for x, y in gene_pairs_sig_tmp]
        
    if 'gene_pairs_ind' in adata.uns:
        gene_pairs_ind_tmp = [tuple(gp.split('_')) for gp in adata.uns['gene_pairs_ind']]
        gene_pairs_ind_tmp = [(x, list(int(val) for val in y.split(' - ')) if ' - ' in y else int(y)) for x, y in gene_pairs_ind_tmp]
        adata.uns['gene_pairs_ind'] = [(list(int(val) for val in x.split(' - ')) if ' - ' in x else int(x), y) for x, y in gene_pairs_ind_tmp]
    
    if 'gene_pairs_sig_ind' in adata.uns:
        gene_pairs_sig_ind_tmp = [tuple(gp.split('_')) for gp in adata.uns['gene_pairs_sig_ind']]
        gene_pairs_sig_ind_tmp = [(x, list(int(val) for val in y.split(' - ')) if ' - ' in y else int(y)) for x, y in gene_pairs_sig_ind_tmp]
        adata.uns['gene_pairs_sig_ind'] = [(list(int(val) for val in x.split(' - ')) if ' - ' in x else int(x), y) for x, y in gene_pairs_sig_ind_tmp]

    def recover_tuple_or_list(g):
        return g.split(' - ') if ' - ' in g else g

    if 'gene_pairs_per_metabolite' in adata.uns:
        for key, subdict in adata.uns['gene_pairs_per_metabolite'].items():
            subdict['gene_pair'] = [
                (recover_tuple_or_list(gp1), recover_tuple_or_list(gp2))
                for gp1, gp2 in subdict['gene_pair']
            ]

    if 'gene_pairs_per_ct_pair' in adata.uns:
        adata.uns['gene_pairs_per_ct_pair'] = {
            tuple(k.split(' - ')): [
                (recover_tuple_or_list(pair[0]), recover_tuple_or_list(pair[1]))
                for pair in v
            ] for k, v in adata.uns['gene_pairs_per_ct_pair'].items()
        }

    if 'gene_pairs_per_ct_pair_ind' in adata.uns:
        adata.uns['gene_pairs_per_ct_pair_ind'] = {
            tuple(k.split(' - ')): v
            for k, v in adata.uns['gene_pairs_per_ct_pair_ind'].items()
        }
    
    if 'gene_pairs_per_ct_pair_sig_ind' in adata.uns:
        adata.uns['gene_pairs_per_ct_pair_sig_ind'] = {
            tuple(k.split(' - ')): v
            for k, v in adata.uns['gene_pairs_per_ct_pair_sig_ind'].items()
        }

    if 'gene_pairs_ind_per_ct_pair' in adata.uns:
        adata.uns['gene_pairs_ind_per_ct_pair'] = {
            tuple(k.split(' - ')): [
                (recover_tuple_or_list(str(pair[0])), recover_tuple_or_list(str(pair[1])))
                for pair in v
            ] for k, v in adata.uns['gene_pairs_ind_per_ct_pair'].items()
        }
    
    if 'gene_pairs_ind_per_ct_pair_sig' in adata.uns:
        adata.uns['gene_pairs_ind_per_ct_pair_sig'] = {
            tuple(k.split(' - ')): [
                (recover_tuple_or_list(str(pair[0])), recover_tuple_or_list(str(pair[1])))
                for pair in v
            ] for k, v in adata.uns['gene_pairs_ind_per_ct_pair_sig'].items()
        }

    def restore_list(x):
        if isinstance(x, str) and ' - ' in x:
            return x.split(' - ')
        return x

    if 'ccc_results' in adata.uns:
        for key in ['cell_com_df_gp', 'cell_com_df_m', 'cell_com_df_gp_sig', 'cell_com_df_m_sig']:
            if key in adata.uns['ccc_results']:
                df = adata.uns['ccc_results'][key]
                adata.uns['ccc_results'][key] = df.applymap(restore_list)
    
    if 'ct_ccc_results' in adata.uns:
        for key in ['cell_com_df_gp', 'cell_com_df_m', 'cell_com_df_gp_sig', 'cell_com_df_m_sig']:
            if key in adata.uns['ct_ccc_results']:
                df = adata.uns['ct_ccc_results'][key]
                adata.uns['ct_ccc_results'][key] = df.applymap(restore_list)

    return


def read_h5ad(
    filename: str,
):
    """
    Load an AnnData object from disk and restore Harreman-specific metadata.

    Parameters
    ----------
    filename : str
        Path to the `.h5ad` file to load.

    Returns
    -------
    AnnData
        The fully restored AnnData object with Harreman metadata recovered.
    """
    
    adata = anndata.read_h5ad(filename)

    if 'genes' in adata.uns:
        adata.uns['genes'] = list(adata.uns['genes'])

    recover_uns_hotspot(adata)
    recover_uns_harreman(adata)

    return adata
