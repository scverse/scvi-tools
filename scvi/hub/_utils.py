import os
from typing import Optional

import anndata


def _get_cell_gene_counts(
    local_dir: str,
    cell_count: Optional[int] = None,
    gene_count: Optional[int] = None,
):
    if os.path.isfile(f"{local_dir}/adata.h5ad"):
        adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
        return adata.n_obs, adata.n_vars
    else:
        if cell_count is None or gene_count is None:
            raise ValueError(
                "No data found on disk. Please provide `cell_count` and `gene_count`."
            )
        else:
            return cell_count, gene_count
