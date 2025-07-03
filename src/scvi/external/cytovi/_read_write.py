import os
from typing import Optional, Union, List, Dict
import anndata as ad
from anndata import AnnData


from pathlib import Path
import numpy as np
import warnings

from scvi.utils import dependencies

@dependencies('readfcs')
def read_fcs(
    path: str,
    return_raw_layer: bool = True,
    include_hidden: bool = False,
    remove_markers: Optional[List[str]] = None,
    sample_name_extractor: Optional[Dict[str, Union[str, int]]] = None,
    verbose: bool = True,
) -> AnnData:
    """Read one or more `.fcs` files and return a concatenated AnnData object.

    This function supports reading a single `.fcs` file or all `.fcs` files in a directory,
    with optional preprocessing including removal of specific markers and extraction of
    sample metadata from file names.

    Parameters
    ----------
    path : str
        Path to a single `.fcs` file or a directory containing `.fcs` files.
    return_raw_layer : bool, optional (default: True)
        Whether to store the untransformed data in the `.layers["raw"]` slot of each AnnData object.
    include_hidden : bool, optional (default: False)
        Whether to include hidden/system files (i.e., those beginning with a period).
    remove_markers : list of str, optional (default: None)
        List of marker names (columns in `.var_names`) to remove from each file.
        Markers not found in the file will be skipped with a warning.
    sample_name_extractor : dict, optional (default: None)
        Dictionary with format `{"split": str, "index": int}`.
        Splits the filename by the provided delimiter and extracts the sample name
        from the specified index position. Example:
        `{"split": "_", "index": 0}` will extract `"sample1"` from `"sample1_batchA.fcs"`.
    verbose : bool, optional (default: True)
        Whether to print progress and warning messages.

    Returns
    -------
    AnnData
        A concatenated AnnData object combining all FCS files, with batch info stored in `.obs`.

    Raises
    ------
    ValueError
        If the provided path is neither a valid `.fcs` file nor a directory.

    Notes
    -----
    - Files with differing marker names will be merged using `join="outer"` (may introduce NaNs).
    - Use the `merge_batches` function to explicitly handle merging of incompatible panels.
    """
    import readfcs

    if path.endswith(".fcs") and os.path.isfile(path):
        fcs_files = [os.path.basename(path)]
        folder = os.path.dirname(path) or "."
    elif os.path.isdir(path):
        folder = path
        fcs_files = [
            f for f in os.listdir(folder)
            if f.endswith('.fcs') and (include_hidden or not f.startswith('.'))
        ]
    else:
        raise ValueError(f"Path '{path}' is neither a valid .fcs file nor a directory.")

    adata_list = []
    var_set_ref = None
    mismatched_files = []

    for i, fcs_file in enumerate(sorted(fcs_files)):
        fcs_path = os.path.join(folder, fcs_file)
        adata = readfcs.read(fcs_path)

        if return_raw_layer:
            adata.layers["raw"] = adata.X.copy()

        if remove_markers:
            present = [m for m in remove_markers if m in adata.var_names]
            missing = [m for m in remove_markers if m not in adata.var_names]
            if missing and verbose:
                print(f"File {fcs_file}: Could not find markers and skipped removal: {missing}")
            if present:
                adata = adata[:, [v for v in adata.var_names if v not in present]]

        adata.obs["filename"] = fcs_file
        adata.obs["sample_id"] = i

        if sample_name_extractor:
            try:
                parts = fcs_file.split(sample_name_extractor["split"])
                sample_name = parts[sample_name_extractor["index"]]
                adata.obs["sample_name"] = sample_name
            except Exception as e:
                if verbose:
                    print(f"Could not extract sample name from '{fcs_file}': {e}")

        adata_list.append(adata)

        var_names = set(adata.var_names)
        if var_set_ref is None:
            var_set_ref = var_names
        elif var_names != var_set_ref:
            mismatched_files.append(fcs_file)

    if mismatched_files and verbose:
        print("Warning: Variable (marker) names differ in the following FCS files:")
        for f in mismatched_files:
            print(f" - {f}")
        print("They will be aligned using `join='outer'`, which may introduce NaNs. "
              "If you want to combine different batches, use the `merge_batches` function.")

    return ad.concat(adata_list, join="outer", label="sample_id", index_unique="-")


@dependencies('fcswrite')
def write_fcs(
    adata: ad.AnnData,
    output_path: str,
    split_by: str | None = None,
    layer: str | None = None,
    verbose: bool = True,
    write_kwargs: Optional[Dict[str, Union[str, int]]] = None,
):
    """
    Export AnnData expression data to one or more FCS files.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, where `adata.X` or `adata.layers[layer]` contains the expression data
        to be written to FCS format.
    output_path : str
        Directory where the FCS files will be written.
    split_by : str, optional
        Column in `adata.obs` to group cells by. If specified, an FCS file will be written for each group.
        If not provided, a single FCS file is written using all cells.
    layer : str, optional
        Layer in `adata.layers` to use as the data source. If `None`, uses `adata.X`.
    verbose : bool, default: True
        If `True`, prints a message for each written file including the output path and data shape.

    Warns
    -----
    UserWarning
        - If non-numeric data is detected in the selected data matrix.
        - If missing values (NaNs) are detected. These rows will be removed before writing.

    Raises
    ------
    ValueError
        If `split_by` is specified but the column is not found in `adata.obs`.

    Notes
    -----
    - This function automatically removes rows with any NaN values before writing to FCS.
    - Data is cast to `float32` before writing.
    - Only `.X` or `.layers[layer]` is written; metadata from `obs` or `var` is not exported to FCS.

    Examples
    --------
    Write a single FCS file:

    >>> write_fcs(adata, "output/")

    Write separate FCS files per condition:

    >>> write_fcs(adata, "output/", split_by="sample_id")

    Use a specific layer:

    >>> write_fcs(adata, "output/", layer="denoised", split_by="batch")

    """
    import fcswrite

    if write_kwargs is None:
        write_kwargs = {}

    output_dir = Path(output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    def get_data(subset, group_name=None):
        data = subset.layers[layer] if layer else subset.X
        if not isinstance(data, np.ndarray):
            data = data.toarray()
        data = data.astype(np.float32)

        # Check for non-numeric data
        if not np.issubdtype(data.dtype, np.number):
            warnings.warn(f"Non-numeric data detected in group '{group_name or 'all'}'. This may cause issues with FCS export.")

        # Check for NaNs before filtering
        n_nans = np.isnan(data).sum()
        if n_nans > 0:
            warnings.warn(f"Group '{group_name or 'all'}' contains {n_nans} NaNs. These rows will be removed before writing.")

        return data

    if split_by:
        if split_by not in adata.obs.columns:
            raise ValueError(f"'{split_by}' not found in adata.obs")
        for group in adata.obs[split_by].unique():
            sub = adata[adata.obs[split_by] == group]
            data = get_data(sub, group)
            data = data[~np.isnan(data).any(axis=1)]
            out_file = output_dir / f"export_{group}.fcs"
            fcswrite.write_fcs(out_file, chn_names=list(sub.var_names), data=data, **write_kwargs)
            if verbose:
                print(f"Wrote: {out_file} | shape={data.shape}")
    else:
        data = get_data(adata)
        data = data[~np.isnan(data).any(axis=1)]
        out_file = output_dir / "export.fcs"
        fcswrite.write_fcs(out_file, chn_names=list(adata.var_names), data=data, **write_kwargs)
        if verbose:
            print(f"Wrote: {out_file} | shape={data.shape}")
