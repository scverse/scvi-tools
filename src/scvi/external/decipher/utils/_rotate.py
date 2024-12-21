import numpy as np
from anndata import AnnData


def _rot(t, u=1):
    if u not in [-1, 1]:
        raise ValueError("u must be -1 or 1")
    return np.array([[np.cos(t), np.sin(t) * u], [-np.sin(t), np.cos(t) * u]])


def rotate_decipher_components(
    adata: AnnData,
    v_obsm_key: str,
    z_obsm_key: str,
    v1_obs_col: str | None = None,
    v1_order: list | None = None,
    v2_obs_col: str | None = None,
    v2_order: list | None = None,
    auto_flip_decipher_z: bool = True,
) -> tuple[np.array, np.array, np.array]:
    """
    Rotate and flip the decipher space v to maximize the correlation of each decipher component.

    The rotation is performed to maximize the correlation of each decipher component with the
    providedcolumns values from `adata.obs` (e.g. pseudo-time, cell state progression, etc.).

    Parameters
    ----------
    adata: sc.AnnData
        The annotated data matrix.
    v_obsm_key: str
        Key in `adata.obsm` to get the decipher v space.
    z_obsm_key: str
        Key in `adata.obsm` to get the decipher z space.
    v1_obs_col: str, optional
        Column name in `adata.obs` to align the first decipher component with.
        If None, only align the second component (or does not align at all if `v2` is also None).
    v1_order: list, optional
        List of values in `adata.obs[v1_col]`. The rotation will attempt to align these values in
        order along the v1 component. Must be provided if `adata.obs[v1_col]` is not numeric.
    v2_obs_col: str, optional
        Column name in `adata.obs` to align the second decipher component with.
        If None, only align the first component (or does not align at all if `v1` is also None).
    v2_order: list, optional
        List of values in `adata.obs[v2_col]`. The rotation will attempt to align these values in
        order along the v2 component. Must be provided if `adata.obs[v2_col]` is not numeric.
    auto_flip_decipher_z: bool, default True
        If True, flip each z to be correlated positively with the components.

    Returns
    -------
    `decipher_v`: ndarray
        The decipher v space after rotation.
    `decipher_z`: ndarray
        The decipher z space after flipping.
    `rotation`: ndarray
        The rotation matrix used to rotate the decipher v space.
    """

    def process_col_obs(v_col, v_order):
        if v_col is not None:
            v_obs = adata.obs[v_col]
            if v_order is not None:
                v_obs = v_obs.astype("category").cat.set_categories(v_order)
                v_obs = v_obs.cat.codes.replace(-1, np.nan)
            elif v_obs.dtype == "category":
                raise ValueError("v_order must be provided if v_obs is a category column")
            v_valid_cells = ~v_obs.isna()
            v_obs = v_obs[v_valid_cells].astype(float)
            return v_obs, v_valid_cells
        return None, None

    v1_obs, v1_valid_cells = process_col_obs(v1_obs_col, v1_order)
    v2_obs, v2_valid_cells = process_col_obs(v2_obs_col, v2_order)

    decipher_v = adata.obsm[v_obsm_key]

    def score_rotation(r):
        rotated_space = decipher_v @ r
        score = 0
        if v1_obs_col is not None:
            score += np.corrcoef(rotated_space[v1_valid_cells, 0], v1_obs)[1, 0]
            score -= np.abs(np.corrcoef(rotated_space[v1_valid_cells, 1], v1_obs)[1, 0])
        if v2_obs_col is not None:
            score += np.corrcoef(rotated_space[v2_valid_cells, 1], v2_obs)[1, 0]
            score -= np.abs(np.corrcoef(rotated_space[v2_valid_cells, 0], v2_obs)[1, 0])
        return score

    if v1_obs_col is None and v2_obs_col is None:
        raise ValueError("At least one of v1_obs_col or v2_obs_col must be provided")

    rotation_scores = []
    for t in np.linspace(0, 2 * np.pi, 100):
        for u in [1, -1]:
            rotation = _rot(t, u)
            rotation_scores.append((score_rotation(rotation), rotation))
    best_rotation = max(rotation_scores, key=lambda x: x[0])[1]

    decipher_v = decipher_v @ best_rotation

    decipher_z = adata.obsm[z_obsm_key]
    if auto_flip_decipher_z:
        # flip each z to be correlated positively with the components
        dim_z = decipher_z.shape[1]
        z_v_corr = np.corrcoef(decipher_z, y=decipher_v, rowvar=False)
        z_sign_correction = np.sign(z_v_corr[:dim_z, dim_z:].sum(axis=1))
        decipher_z = decipher_z * z_sign_correction

    return decipher_v, decipher_z, best_rotation
