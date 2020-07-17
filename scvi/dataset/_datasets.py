from typing import Optional
import anndata

from ._synthetic import _generate_synthetic


def synthetic_iid(
    batch_size: Optional[int] = 200,
    n_genes: Optional[int] = 100,
    n_proteins: Optional[int] = 100,
    n_batches: Optional[int] = 2,
    n_labels: Optional[int] = 3,
) -> anndata.AnnData:
    """Synthetic dataset with ZINB distributed RNA and NB distributed protein

    Each value is independently and identically distributed.

    Parameters
    ----------
    batch_size
        Number of cells per batch
    n_genes
        Number of genes
    n_proteins
        Number of proteins
    n_batches
        Number of batches
    n_labels
        Number of cell types

    Returns
    -------
    `AnnData` object

    """

    return _generate_synthetic(
        batch_size=batch_size,
        n_genes=n_genes,
        n_proteins=n_proteins,
        n_batches=n_batches,
        n_labels=n_labels,
    )
