import anndata
import pooch

from scvi.external.harreman._data import harreman_data_hash, harreman_data_url

DATASETS_DIR = "datasets"


def load_visium_mouse_colon_dataset(
    sample: str | None = None,
    save_path: str = "data/",
) -> anndata.AnnData:
    """
    Load the mouse colon 10x Visium dataset (Parigi et al.).

    Parameters
    ----------
    sample : str, optional
        Sample to load. One of "d0", "d14". If None, loads the full unrolled dataset.
    save_path : str, optional
        Path to save the downloaded dataset.

    Returns
    -------
    adata : AnnData
        The loaded 10x Visium dataset.
    """
    dataset_prefix = "Parigi_et_al_mouse_colon"

    samples_info = {
        "d0": f"{dataset_prefix}_d0.h5ad",
        "d14": f"{dataset_prefix}_d14.h5ad",
    }

    if sample is not None:
        if sample not in samples_info:
            raise ValueError(f'"sample" needs to be one of: {list(samples_info.keys())}')
        fname = samples_info[sample]
    else:
        fname = f"{dataset_prefix}_unrolled.h5ad"

    path = pooch.retrieve(
        url=harreman_data_url(DATASETS_DIR, fname),
        known_hash=harreman_data_hash(DATASETS_DIR, fname),
        fname=fname,
        path=save_path,
        progressbar=True,
    )

    return anndata.read_h5ad(path)


def load_slide_seq_human_lung_dataset(
    sample: str | None = None,
    save_path: str = "data/",
) -> anndata.AnnData:
    """
    Load the human lung Slide-seq dataset (Liu et al.).

    Parameters
    ----------
    sample : str, optional
        Sample to load. One of "Puck_200727_08", "Puck_200727_09", "Puck_200727_10",
        "Puck_220408_13", "Puck_220408_14", "Puck_220408_15", "Puck_220408_20".
        If None, loads the full combined dataset.
    save_path : str, optional
        Path to save the downloaded dataset.

    Returns
    -------
    adata : AnnData
        The loaded Slide-seq dataset.
    """
    dataset_prefix = "Liu_et_al_human_lung"

    samples_info = {
        "Puck_200727_08": f"{dataset_prefix}_Puck_200727_08.h5ad",
        "Puck_200727_09": f"{dataset_prefix}_Puck_200727_09.h5ad",
        "Puck_200727_10": f"{dataset_prefix}_Puck_200727_10.h5ad",
        "Puck_220408_13": f"{dataset_prefix}_Puck_220408_13.h5ad",
        "Puck_220408_14": f"{dataset_prefix}_Puck_220408_14.h5ad",
        "Puck_220408_15": f"{dataset_prefix}_Puck_220408_15.h5ad",
        "Puck_220408_20": f"{dataset_prefix}_Puck_220408_20.h5ad",
    }

    if sample is not None:
        if sample not in samples_info:
            raise ValueError(f'"sample" needs to be one of: {list(samples_info.keys())}')
        fname = samples_info[sample]
    else:
        fname = f"{dataset_prefix}.h5ad"

    path = pooch.retrieve(
        url=harreman_data_url(DATASETS_DIR, fname),
        known_hash=harreman_data_hash(DATASETS_DIR, fname),
        fname=fname,
        path=save_path,
        progressbar=True,
    )

    return anndata.read_h5ad(path)
