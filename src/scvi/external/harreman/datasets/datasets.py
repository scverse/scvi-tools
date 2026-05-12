
import os
import tempfile
from typing import Optional
from pathlib import Path
import scanpy as sc


temp_dir_obj = tempfile.TemporaryDirectory()


def load_visium_mouse_colon_dataset(
    sample: Optional[str] = None,
) -> "sc.AnnData":
    """
    Load the mouse colon 10x Visium dataset.

    Returns
    -------
    adata : AnnData
        The loaded 10x Visium dataset.
    """
    
    dataset_prefix = 'Parigi_et_al_mouse_colon'
    
    samples_path_dict = {
        'd0': 'https://figshare.com/ndownloader/files/59325113',
        'd14': 'https://figshare.com/ndownloader/files/59325116',
    }
    
    if sample:
        if sample not in samples_path_dict.keys():
            raise ValueError(f'"sample" needs to be one of: {list(samples_path_dict.keys())}')
        else:
            adata_path = os.path.join(temp_dir_obj.name, f"{dataset_prefix}_{sample}.h5ad")
            backup_url = samples_path_dict[sample]
    else:
        adata_path = os.path.join(temp_dir_obj.name, f"{dataset_prefix}_unrolled.h5ad")
        backup_url = 'https://figshare.com/ndownloader/files/59325119'
    
    adata = sc.read(adata_path, backup_url=backup_url)

    return adata


def load_slide_seq_human_lung_dataset(
    sample: Optional[str] = None,
) -> "sc.AnnData":
    """
    Load the human lung Slide-seq dataset.

    Returns
    -------
    adata : AnnData
        The loaded Slide-seq dataset.
    """
    
    dataset_prefix = 'Liu_et_al_human_lung'
    
    samples_path_dict = {
        'Puck_200727_08': 'https://figshare.com/ndownloader/files/59325098',
        'Puck_200727_09': 'https://figshare.com/ndownloader/files/59325092',
        'Puck_200727_10': 'https://figshare.com/ndownloader/files/59325095',
        'Puck_220408_13': 'https://figshare.com/ndownloader/files/59325101',
        'Puck_220408_14': 'https://figshare.com/ndownloader/files/59325104',
        'Puck_220408_15': 'https://figshare.com/ndownloader/files/59325107',
        'Puck_220408_20': 'https://figshare.com/ndownloader/files/59325110',
    }

    if sample:
        if sample not in samples_path_dict.keys():
            raise ValueError(f'"sample" needs to be one of: {list(samples_path_dict.keys())}')
        else:
            adata_path = os.path.join(temp_dir_obj.name, f"{dataset_prefix}_{sample}.h5ad")
            backup_url = samples_path_dict[sample]
    else:
        adata_path = os.path.join(temp_dir_obj.name, f"{dataset_prefix}.h5ad")
        backup_url = 'https://figshare.com/ndownloader/files/59325125'
    
    adata = sc.read(adata_path, backup_url=backup_url)

    return adata
