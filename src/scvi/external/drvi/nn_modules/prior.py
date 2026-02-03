from __future__ import annotations

import torch
from torch import nn
from torch.distributions import Normal, kl_divergence


class Prior(nn.Module):
    """Abstract base class for prior distributions in variational autoencoders.

    This class defines the interface for different prior distributions that can be
    used in variational autoencoders. The prior distribution is used to regularize
    the latent space and compute the KL divergence with the posterior.

    Notes
    -----
    Subclasses must implement:
    - `kl`: Method to compute KL divergence with the posterior distribution
    """

    def __init__(self) -> None:
        super().__init__()

    def kl(self, qz: Normal) -> torch.Tensor:
        """Compute KL divergence between posterior and prior.

        Parameters
        ----------
        qz
            The posterior distribution (usually Normal).

        Returns
        -------
        torch.Tensor
            KL divergence between posterior and prior.

        Raises
        ------
        NotImplementedError
            This method must be implemented by subclasses.
        """
        raise NotImplementedError()


class StandardPrior(Prior):
    """Standard normal prior distribution.

    This is the most commonly used prior in variational autoencoders.
    It assumes a standard normal distribution N(0, I) for the latent variables.

    Examples
    --------
    >>> import torch
    >>> from torch.distributions import Normal
    >>> # Create standard prior
    >>> prior = StandardPrior()
    >>> # Create a posterior distribution
    >>> qz = Normal(torch.randn(10, 5), torch.ones(10, 5))
    >>> # Compute KL divergence
    >>> kl = prior.kl(qz)
    >>> print(kl.shape)  # torch.Size([10, 5])
    """

    def __init__(self) -> None:
        super().__init__()

    def kl(self, qz: Normal) -> torch.Tensor:
        """Compute KL divergence with standard normal prior.

        Parameters
        ----------
        qz
            The posterior distribution (must be Normal).

        Returns
        -------
        torch.Tensor
            KL divergence between qz and N(0, I).

        Notes
        -----
        For Normal distributions, the KL divergence has a closed-form solution:
        KL(q||p) = 0.5 * (μ² + σ² - log(σ²) - 1)
        where μ and σ are the mean and standard deviation of q.
        """
        # 1 x N
        assert isinstance(qz, Normal)
        return kl_divergence(qz, Normal(torch.zeros_like(qz.mean), torch.ones_like(qz.mean)))
