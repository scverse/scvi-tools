import os
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
from anndata import AnnData

from scvi import settings
from scvi.utils import dependencies


@dependencies("readfcs")
def read_fcs(
    path: str,
    return_raw_layer: bool = True,
    include_hidden: bool = False,
    remove_markers: list[str] | None = None,
    sample_name_extractor: dict[str, str | int] | None = None,
    verbose: bool = True,
) -> AnnData:
    """Read one or more `.fcs` files and return a concatenated AnnData object.

    This function supports reading a single `.fcs` file or all `.fcs` files in a directory,
    with optional preprocessing including removal of specific markers and extraction of
    sample metadata from file names.

    Parameters
    ----------
    path
        Path to a single `.fcs` file or a directory containing `.fcs` files.
    return_raw_layer
        Whether to store the untransformed data in the `.layers["raw"]` slot of each AnnData
        object. Default: True.
    include_hidden
        Whether to include hidden/system files (i.e., those beginning with a period).
        Default: False.
    remove_markers
        List of marker names (columns in `.var_names`) to remove from each file.
        Markers not found in the file will be skipped with a warning. Default: None.
    sample_name_extractor
        Dictionary with format `{"split": str, "index": int}`.
        Splits the filename by the provided delimiter and extracts the sample name
        from the specified index position. Example:
        `{"split": "_", "index": 0}` will extract `"sample1"` from `"sample1_batchA.fcs"`.
        Default: None.
    verbose
        Whether to print progress and warning messages. Default: True.

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
            f
            for f in os.listdir(folder)
            if f.endswith(".fcs") and (include_hidden or not f.startswith("."))
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
            except (IndexError, KeyError, ValueError) as e:
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
        print(
            "They will be aligned using `join='outer'`, which may introduce NaNs. "
            "If you want to combine different batches, use the `merge_batches` function."
        )

    return ad.concat(adata_list, join="outer", label="sample_id", index_unique="-")


@dependencies("fcswrite")
def write_fcs(
    adata: ad.AnnData,
    output_path: str | None = None,
    split_by: str | None = None,
    layer: str | None = None,
    prefix: str = "export",
    verbose: bool = True,
    write_kwargs: dict[str, str | int] | None = None,
):
    """
    Export AnnData expression data to one or more FCS files.

    Parameters
    ----------
    adata
        Annotated data matrix, where `adata.X` or `adata.layers[layer]` contains the expression
        data to be written to FCS format.
    output_path
        Directory where the FCS files will be written. If `None`, uses the current working
        directory.
    split_by
        Column in `adata.obs` to group cells by. If specified, an FCS file will be written for
        each group. If not provided, a single FCS file is written using all cells.
    layer
        Layer in `adata.layers` to use as the data source. If `None`, uses `adata.X`.
    prefix
        Prefix for the output FCS file names. Default: "export".
    verbose
        If `True`, prints a message for each written file including the output path and data shape.
        Default: True.
    write_kwargs
        Additional keyword arguments to pass to `fcswrite.write_fcs`. This can include
        parameters like `text_kw_pr` for custom text annotations in the FCS file.

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
    - Only `.X` or `.layers[layer]` is written; metadata from `obs` or `var` is not exported
      to FCS.

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

    output_dir = Path(output_path) if output_path else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)

    chn_names = list(adata.var_names)
    text = {f"$P{i + 1}N": chn for i, chn in enumerate(chn_names)} | {
        f"$P{i + 1}S": chn for i, chn in enumerate(chn_names)
    }

    def get_data(subset, group_name=None):
        data = subset.layers[layer] if layer else subset.X
        if not isinstance(data, np.ndarray):
            data = data.toarray()
        data = data.astype(np.float32)

        # Check for non-numeric data
        if not np.issubdtype(data.dtype, np.number):
            warnings.warn(
                f"Non-numeric data detected in group '{group_name or 'all'}'. This may cause "
                "issues with FCS export.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        # Check for NaNs before filtering
        n_nans = np.isnan(data).sum()
        if n_nans > 0:
            warnings.warn(
                f"Group '{group_name or 'all'}' contains {n_nans} NaNs. These rows will be "
                "removed before writing.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        return data

    if split_by:
        if split_by not in adata.obs.columns:
            raise ValueError(f"'{split_by}' not found in adata.obs")
        for group in adata.obs[split_by].unique():
            sub = adata[adata.obs[split_by] == group]
            data = get_data(sub, group)
            data = data[~np.isnan(data).any(axis=1)]
            out_file = output_dir / f"{prefix}_{group}.fcs"
            fcswrite.write_fcs(
                out_file,
                chn_names=chn_names,
                data=data,
                text_kw_pr=text,
                compat_percent=False,
                **write_kwargs,
            )
            if verbose:
                print(f"Wrote: {out_file} | shape={data.shape}")
    else:
        data = get_data(adata)
        data = data[~np.isnan(data).any(axis=1)]
        out_file = output_dir / f"{prefix}.fcs"
        fcswrite.write_fcs(
            out_file,
            chn_names=chn_names,
            data=data,
            text_kw_pr=text,
            compat_percent=False,
            **write_kwargs,
        )
        if verbose:
            print(f"Wrote: {out_file} | shape={data.shape}")
