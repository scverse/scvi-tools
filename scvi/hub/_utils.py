import os
from typing import Optional

import anndata

from scvi.data._utils import _is_latent


# TODO remove if unused
def _get_counts_and_latent_info(
    local_dir: str,
    cell_count: Optional[int] = None,
    gene_count: Optional[int] = None,
    is_latent: Optional[bool] = None,
):
    if os.path.isfile(f"{local_dir}/adata.h5ad"):
        adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
        return adata.n_obs, adata.n_vars, _is_latent(adata)
    else:
        if cell_count is None or gene_count is None:
            raise ValueError(
                "No data found on disk. Please provide `cell_count` and `gene_count`."
            )
        else:
            return cell_count, gene_count, is_latent
