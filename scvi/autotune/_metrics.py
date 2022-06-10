<<<<<<< HEAD
=======
from typing import Optional

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples, silhouette_score

from scvi.model.base import BaseModelClass


def silhouette_metric_labels_batch(
    model: BaseModelClass,
    labels_key: str,
    batch_key: str,
    sample_size: int = 1000,
<<<<<<< HEAD
) -> float:
    """
    Batch- and label-wise silhouette.

=======
    unknown_label: Optional[str] = None,
) -> float:
    """
    Batch- and label-wise silhouette.
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    Parameters
    ----------
    model
        A scvi model.
    labels_key
        The key of the labels.
    batch_key
        The key of the batch.
    sample_size
        Sample size for silhouette. Randomly subsets the data
<<<<<<< HEAD

=======
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    Returns
    -------
    Sum of asw for batch and for labels. Scores are scaled such that
    2 is the best score and 0 is the worst.
<<<<<<< HEAD

    Notes
    -----
    This function is influence by the following code:

=======
    Notes
    -----
    This function is influence by the following code:
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    https://github.com/theislab/scib/blob/main/scib/metrics/silhouette.py
    """
    model.is_trained_ = True
    latent = model.get_latent_representation()
    model.is_trained_ = False
    adata = model.adata

<<<<<<< HEAD
=======
    if unknown_label:
        adata = adata[adata.obs[labels_key] != unknown_label].copy()

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    # bio conservation
    asw_labels = silhouette_score(
        latent,
        adata.obs[labels_key],
        sample_size=sample_size,
    )
    # normalize into 0-1 range
    asw_labels = (asw_labels + 1) / 2

    sil_all = pd.DataFrame(columns=["group", "silhouette_score"])

    for group in adata.obs[labels_key].unique():
        mask = np.asarray(adata.obs[labels_key] == group)
        adata_group = adata[mask]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil_per_group = silhouette_samples(latent[mask], adata_group.obs[batch_key])

        # take only absolute value
        sil_per_group = np.abs(sil_per_group)

        # scale it
        sil_per_group = 1 - sil_per_group

        sil_all = sil_all.append(
            pd.DataFrame(
                {
                    "group": [group] * len(sil_per_group),
                    "silhouette_score": sil_per_group,
                }
            )
        )

    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby("group").mean()
    asw_batch = sil_means["silhouette_score"].mean()

    return asw_batch + asw_labels
