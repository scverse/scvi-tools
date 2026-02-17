import numpy as np
import scipy.spatial


def compute_foscttm(
    latents: dict[str, np.ndarray],
    indices: list[np.ndarray] | None = None,
    downsample: bool = False,
    n_obs: int = 10000,
) -> dict[str, float]:
    """Compute FOSCTTM (Fraction of Samples Closer Than True Match) on latent embeddings.

    For paired multi-modal data, measures how well the model aligns
    corresponding cells in the latent space. Lower values are better
    (0 = perfect alignment). Assumes cells are aligned: mod1 cell i 
    corresponds to mod2 cell i.

    Parameters
    ----------
    latents
        Dictionary mapping modality names to latent arrays.
    indices
        List of two arrays: indices for [mod1, mod2].
        Usually obtained from model.validation_indices after training.
    downsample
        Whether to downsample when the number of observations is large.
    n_obs
        Number of observations to keep if downsample is True.
    """
    # Extract latent representations for the two modalities
    latent_mod1 = latents[list(latents.keys())[0]]
    latent_mod2 = latents[list(latents.keys())[1]]

    # Extract indices for the two modalities if provided
    if indices is not None:
        latent_mod1 = latent_mod1[indices[0]]
        latent_mod2 = latent_mod2[indices[1]]

    # Validate shapes match
    if latent_mod1.shape[1] != latent_mod2.shape[1]:
        raise ValueError("Shapes do not match!")

    n_cells = latent_mod1.shape[0]

     # Downsample if requested and number of cells exceeds n_obs
    if n_cells > n_obs and downsample:
        np.random.seed(0)
        sample_indices = np.random.choice(n_cells, size=n_obs, replace=False)
        latent_mod1 = latent_mod1[sample_indices]
        latent_mod2 = latent_mod2[sample_indices]
        n_cells = n_obs

    # Compute pairwise distances
    distances = scipy.spatial.distance_matrix(latent_mod1, latent_mod2)

    # Compute FOSCTTM
    foscttm_mod1_to_mod2 = (distances < np.expand_dims(np.diag(distances), axis=1)).mean(axis=1)
    foscttm_mod2_to_mod1 = (distances < np.expand_dims(np.diag(distances), axis=0)).mean(axis=0)

    # Aggregate metrics
    mean_mod1_to_mod2 = foscttm_mod1_to_mod2.mean()
    mean_mod2_to_mod1 = foscttm_mod2_to_mod1.mean()
    foscttm_metrics = {
        "foscttm/mod1_to_mod2": float(mean_mod1_to_mod2),
        "foscttm/mod2_to_mod1": float(mean_mod2_to_mod1),
        "foscttm/mean": float((mean_mod1_to_mod2 + mean_mod2_to_mod1) / 2),
    }

    return foscttm_metrics
