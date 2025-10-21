import numpy as np
from rich import print
from scipy.sparse import csr_matrix

from scvi.utils import dependencies


def adjusted_nearest_neighbors(
    cell_samples: np.array,
    cell_coordinates: np.array,
    cell_labels: np.array,
    radius: int | None = None,
    k_nn: int | None = None,
    return_sparse: bool = True,
) -> np.ndarray | csr_matrix:
    """Compute the adjacency matrix for the neighborhood of each sample."

    Parameters
    ----------
    cell_samples
        Array of sample id for each cell.
    cell_coordinates
        Array of cell coordinates.
    cell_labels
        Array of cell labels.
    radius
        Radius for the nearest neighbors search.
    k_nn
        Number of nearest neighbors to consider.
    return_sparse
        Whether to return a sparse matrix.

    Returns
    -------
    np.ndarray | csr_matrix
        Adjacency matrix for the neighborhood of each sample,
        without connections to cells of the same type.
    """
    from scipy.sparse import block_diag
    from sklearn.neighbors import NearestNeighbors

    adjacency_matrices = []

    sample_names = np.unique(cell_samples)
    print(f"Using {len(sample_names)} samples")

    if (radius is None) and (k_nn is None):
        raise ValueError("Either radius or k_nn must be provided")

    for sample in sample_names:
        mask = np.squeeze(cell_samples == sample, axis=1)  # n_cells
        sample_coord = cell_coordinates[mask]  # n_cells_sample_i x 2
        sample_cell_types = np.squeeze(cell_labels[mask], axis=1)  # n_cells_sample_i

        if radius is not None:
            nn = NearestNeighbors(radius=radius)
            nn.fit(sample_coord)
            A = nn.radius_neighbors_graph(sample_coord)
        elif k_nn is not None:
            nn = NearestNeighbors(n_neighbors=k_nn + 1)
            nn.fit(sample_coord)
            A = nn.kneighbors_graph(sample_coord)
        else:
            raise ValueError("Either radius or k_nn must be provided.")

        # Find rows and columns of non-zero entries
        row_indices, col_indices = A.nonzero()

        # Create a sparse mask of entries to zero out
        # Only zero out entries where row and column have the same label
        mask_matrix = np.where(sample_cell_types[row_indices] == sample_cell_types[col_indices])[0]

        A_adjusted = A.copy()
        # Zero out these specific data entries
        A_adjusted.data[mask_matrix] = 0

        A_adjusted.eliminate_zeros()

        adjacency_matrices.append(A_adjusted.astype(bool, copy=False))

    adjacency_matrix = block_diag(adjacency_matrices, format="csr")

    row_counts = np.diff(adjacency_matrix.indptr)
    # print mean and std of number of neighbors, round to 2 decimals:
    print(f"Mean number of neighbors: {np.mean(row_counts):.1f} ± {np.std(row_counts):.1f}")

    if radius is None:
        assert np.mean(row_counts) <= k_nn, "Mean number of neighbors is greater than k_nn"

    if return_sparse:
        return adjacency_matrix

    return adjacency_matrix.toarray()


def _get_nonzero_indices_from_rows(csr_matrix, row_idx):
    return np.unique(csr_matrix[row_idx].indices)


@dependencies("matplotlib")
def get_connectivity_distribution(csr_matrix):
    import matplotlib.pyplot as plt

    # Get the number of non-zero entries per row
    row_counts = np.diff(csr_matrix.indptr)

    # Create the histogram
    fig, ax = plt.subplots()
    ax.hist(row_counts, bins=np.max(row_counts) + 1, density=True)
    ax.set_xlabel("Number of non-zero entries")
    ax.set_ylabel("Number of rows")
    ax.set_title("Histogram of non-zero entries per row")
    plt.show()
