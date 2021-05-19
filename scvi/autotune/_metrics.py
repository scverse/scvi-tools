import scanpy as sc
import torch
from sklearn.metrics import silhouette_score


def silhouette_metric(model):
    model.is_trained_ = True
    latent = model.get_latent_representation()
    model.is_trained_ = False
    adata = model.adata
    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.5)
    return silhouette_score(
        adata.obsm["X_scVI"],
        adata.obs["leiden_scVI"],
    )


def _compute_fast_mmd(z1, z2):
    """
    Computes a fast approximation of the MMD. See `_compute_mmd`.
    Parameters
    ----------
    z1
        First tensor
    z2
        Second tensor
    Returns
    -------
    Tensor with one item containing an approximation of the MMD of ``z1`` and ``z2``.
    """
    z1_size = z1.size(0)
    z2_size = z2.size(0)

    # z1_size and z2_size must match, otherwise pick their min
    batch_size = min(z1_size, z2_size)
    print(batch_size)

    # Drop a sample if batch_size is not even
    if batch_size % 2 != 0:
        batch_size -= 1

    z1 = z1[:batch_size, :]
    z2 = z2[:batch_size, :]

    if z1.size(1) != z2.size(1):
        raise ValueError(
            "z1 and z2 must be defined on the same space, "
            "but input was: z1_d={} while "
            "z2_d={}.".format(z1.size(1), z2.size(1))
        )

    # zi_even contains every other sample in zi starting at 0
    z1_even = z1[: batch_size - 1 : 2, :]
    z1_odd = z1[1:batch_size:2, :]
    z2_even = z2[: batch_size - 1 : 2, :]
    z2_odd = z2[1:batch_size:2, :]

    # Compute the kernels
    z1_even_odd_kernels = torch.exp(-((z1_even - z1_odd).pow(2).sum(1)))
    z2_even_odd_kernels = torch.exp(-((z2_even - z2_odd).pow(2).sum(1)))
    z1_even_z2_odd_kernels = torch.exp(-((z1_even - z2_odd).pow(2).sum(1)))
    z1_odd_z2_even_kernels = torch.exp(-((z1_odd - z2_even).pow(2).sum(1)))

    all_kernels = (
        z1_even_odd_kernels
        + z2_even_odd_kernels
        - z1_even_z2_odd_kernels
        - z1_odd_z2_even_kernels
    )
    mmd = all_kernels.mean()
    return mmd


def _compute_mmd_loss(z: torch.Tensor, batch_indices: torch.Tensor) -> torch.Tensor:
    """
    Computes the overall MMD loss associated with this set of samples. The overall MMD is the sum
    of batch-wise MMD's, i.e. the MMD associated with the samples from ``z`` for each pair of sequential
    batches in batch_indices.
    Parameters
    ----------
    z
        Set of samples to compute the overall MMD loss on
    batch_indices
        Batch indices corresponding to each sample in ``z``. Same length as ``z``.
    Returns
    -------
    Tensor with one item containing the MMD loss for the samples in ``z``
    """
    mmd_loss = 0.0
    batches = torch.unique(batch_indices)
    for b0, b1 in zip(batches, batches[1:]):
        z0 = z[(batch_indices == b0).reshape(-1)]
        z1 = z[(batch_indices == b1).reshape(-1)]

        mmd_loss += _compute_fast_mmd(z0, z1)

    return mmd_loss


def mmd(model):
    return _compute_mmd_loss(
        torch.tensor(model.adata.obsm["X_scVI"]),
        torch.tensor(model.adata.obs._scvi_batch),
    )
