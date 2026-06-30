"""Utility functions for joint embedding loss."""

from __future__ import annotations

import torch


def variance_loss(z: torch.Tensor, target_std: float = 1.0) -> torch.Tensor:
    """Compute variance regularization loss (VICReg-style).

    Penalizes dimensions with standard deviation below target_std.
    This prevents collapse where all samples map to similar embeddings.

    Parameters
    ----------
    z
        Embedding tensor of shape (batch_size, n_latent).
    target_std
        Target standard deviation for each dimension. Default is 1.0.

    Returns
    -------
    torch.Tensor
        Scalar loss value (mean of hinge losses across dimensions).
    """
    # Use the biased estimator so the loss is well-defined (0, not NaN) for
    # single-sample minibatches, which the default data splitter can emit.
    std = z.std(dim=0, unbiased=False)
    return torch.relu(target_std - std).mean()


def binomial_split(
    x: torch.Tensor, p: float | torch.Tensor = 0.5
) -> tuple[torch.Tensor, torch.Tensor]:
    """Split counts via binomial thinning.

    Parameters
    ----------
    x
        Input count tensor of shape (batch_size, n_features).
    p
        Probability for binomial split. Can be a scalar or a tensor of shape
        (batch_size, 1) for per-cell probabilities. Default is 0.5.

    Returns
    -------
    tuple
        Two tensors (x1, x2) where x1 ~ Binomial(x, p) and x2 = x - x1.
        Guarantees x1 + x2 = x and both are non-negative.

    Notes
    -----
    Binomial thinning requires non-negative integer counts. Inputs are rounded
    and clamped to ``>= 0`` so that data stored as integer-valued floats (the
    usual single-cell convention) is handled, and non-count inputs do not raise
    mid-training.
    """
    counts = x.round().clamp(min=0)
    x1 = torch.distributions.Binomial(total_count=counts, probs=p).sample()
    x2 = counts - x1
    return x1, x2


def sample_thinning_probs(x: torch.Tensor, min_library_size: float = 10.0) -> torch.Tensor:
    """Sample per-cell thinning probabilities for realistic library size variation.

    Samples target library sizes from a log-uniform distribution between
    min_library_size and the observed library size for each cell. This
    matches the variation typically observed in real single-cell data.

    Parameters
    ----------
    x
        Input count tensor of shape (batch_size, n_features).
    min_library_size
        Minimum target library size. Default is 10.

    Returns
    -------
    torch.Tensor
        Per-cell thinning probabilities of shape (batch_size, 1).
    """
    # Compute observed library sizes per cell
    lib_size = x.sum(dim=1, keepdim=True)  # (batch_size, 1)

    # Sample target library size from LogUniform(min_library_size, lib_size)
    # log(L_target) ~ Uniform(log(min_lib), log(lib_size))
    log_min = torch.log(torch.tensor(min_library_size, device=x.device, dtype=x.dtype))
    log_lib = torch.log(lib_size.clamp(min=min_library_size + 1e-6))

    # Sample uniformly in log space
    log_target = log_min + torch.rand_like(lib_size) * (log_lib - log_min)
    target_lib_size = torch.exp(log_target)

    # Compute thinning probability: p = target_lib_size / lib_size
    # Clamp to [0, 1] for safety (handles edge cases where lib_size <= min_library_size)
    p = (target_lib_size / lib_size.clamp(min=1e-6)).clamp(0.0, 1.0)

    return p


def cross_correlation_loss(
    z1: torch.Tensor,
    z2: torch.Tensor,
    lambda_off_diag: float = 0.01,
    return_components: bool = False,
) -> torch.Tensor | dict[str, torch.Tensor]:
    """Compute cross-correlation objective (CCO) loss.

    This loss has two components:
    - Diagonal (invariance): encourages C_ii = 1, meaning corresponding dimensions
      of the two views should be correlated. This is the self-supervision signal.
    - Off-diagonal (redundancy reduction): encourages C_ij = 0 for i != j, meaning
      different latent dimensions should be decorrelated. This prevents mode collapse.

    Parameters
    ----------
    z1
        First embedding tensor of shape (batch_size, n_latent).
    z2
        Second embedding tensor of shape (batch_size, n_latent).
    lambda_off_diag
        Weight for off-diagonal penalty. Default is 0.01.
    return_components
        If True, return a dict with individual components. Default is False.

    Returns
    -------
    torch.Tensor or dict
        If return_components is False, returns scalar loss value.
        If True, returns dict with keys: 'total', 'invariance', 'redundancy'.
    """
    # Normalize (z-score across batch dimension). The biased std keeps this
    # finite for single-sample minibatches and makes C a proper correlation matrix.
    z1_norm = (z1 - z1.mean(0)) / (z1.std(0, unbiased=False) + 1e-8)
    z2_norm = (z2 - z2.mean(0)) / (z2.std(0, unbiased=False) + 1e-8)

    # Cross-correlation matrix: C[i,j] = corr(z1[:,i], z2[:,j])
    batch_size = z1.shape[0]
    C = z1_norm.T @ z2_norm / batch_size

    # Diagonal loss (invariance): (1 - C_ii)^2 - self-supervision signal
    invariance_loss = ((1 - C.diagonal()) ** 2).sum()

    # Off-diagonal loss (redundancy reduction): sum(C_ij^2 for i != j) - anti-collapse
    # Clone C before modifying diagonal to avoid in-place operation issues
    C_off_diag = C.clone()
    C_off_diag.fill_diagonal_(0)
    redundancy_loss = (C_off_diag**2).sum()

    total_loss = invariance_loss + lambda_off_diag * redundancy_loss

    if return_components:
        return {
            "total": total_loss,
            "invariance": invariance_loss,
            "redundancy": redundancy_loss,
        }
    return total_loss
