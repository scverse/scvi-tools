from __future__ import annotations

from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.distributions import Normal, kl_divergence

from scvi.external.drvi.module._constants import MODULE_KEYS

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

# Standard, VaMP, GMM from Karin's CSI repo


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


class VampPrior(Prior):
    """VampPrior (Variational Mixture of Posteriors) prior distribution.

    This prior uses a mixture of posteriors as the prior distribution,
    which can capture more complex structure in the data than a standard
    normal prior.

    This is adapted from https://github.com/jmtomczak/intro_dgm/main/vaes/vae_priors_example.ipynb

    Parameters
    ----------
    n_components
        Number of mixture components.
    encoder
        Encoder network to compute posterior parameters for pseudo-inputs.
    model_input
        Dictionary containing input data for pseudo-inputs.
    trainable_keys
        Keys in model_input that should be trainable parameters.
    fixed_keys
        Keys in model_input that should be fixed parameters.
    input_type
        Type of input format expected by the encoder.
    preparation_function
        Function to prepare inputs for the encoder.

    Notes
    -----
    VampPrior learns a set of pseudo-inputs and uses the posterior
    distributions of these pseudo-inputs as mixture components for
    the prior. This allows the prior to adapt to the data distribution.

    The prior is computed as:
    p(z) = Σ_k w_k * q(z|u_k)
    where u_k are the pseudo-inputs and w_k are mixing weights.

    Examples
    --------
    >>> import torch
    >>> from torch import nn
    >>> # Create a simple encoder
    >>> class SimpleEncoder(nn.Module):
    ...     def __init__(self):
    ...         super().__init__()
    ...         self.fc = nn.Linear(10, 4)
    ...
    ...     def forward(self, x):
    ...         return self.fc(x), torch.ones_like(self.fc(x))
    >>> # Create model input
    >>> model_input = {"x": torch.randn(5, 10)}
    >>> # Create VampPrior
    >>> prior = VampPrior(n_components=5, encoder=SimpleEncoder(), model_input=model_input, trainable_keys=("x",))
    """

    # Adapted from https://github.com/jmtomczak/intro_dgm/main/vaes/vae_priors_example.ipynb
    # K - components, I - inputs, L - latent, N - samples
    def __init__(
        self,
        n_components: int,
        encoder: nn.Module,
        model_input: dict[str, Any],
        trainable_keys: tuple[str, ...] = ("x",),
        fixed_keys: tuple[str, ...] = (),
        input_type: str = "scvi",
        preparation_function: Callable | None = None,
    ) -> None:
        super().__init__()

        self.encoder = encoder
        self.input_type = input_type
        self.preparation_function = preparation_function

        # pseudo inputs
        pi_aux_data = {}
        pi_tensor_data = {}
        for key in fixed_keys:
            if isinstance(model_input[key], torch.Tensor):
                pi_tensor_data[key] = torch.nn.Parameter(model_input[key], requires_grad=False)
            else:
                pi_aux_data[key] = model_input[key]
        for key in trainable_keys:
            pi_tensor_data[key] = torch.nn.Parameter(model_input[key], requires_grad=True)
            assert pi_tensor_data[key].shape[0] == n_components
        self.pi_aux_data = pi_aux_data
        self.pi_tensor_data = nn.ParameterDict(pi_tensor_data)

        # mixing weights
        self.w = torch.nn.Parameter(torch.zeros(n_components, 1, 1))  # K x 1 x 1

    def get_params(self) -> tuple[torch.Tensor, torch.Tensor]:
        """Get the parameters of the mixture components.

        Returns
        -------
        tuple
            (means, variances) of the mixture components.

        Notes
        -----
        This method computes the posterior parameters for each pseudo-input
        by passing them through the encoder network.
        """
        # u->encoder->mean, var
        original_mode = self.encoder.training
        self.encoder.train(False)
        if self.input_type == "scfemb":
            z = self.encoder({**self.pi_aux_data, **self.pi_tensor_data})
            output = z[MODULE_KEYS.QZM_KEY], z[MODULE_KEYS.QZV_KEY]
        elif self.input_type == "scvi":
            if self.preparation_function is None:
                raise ValueError("preparation_function must be provided for scvi input type")
            x, args, kwargs = self.preparation_function({**self.pi_aux_data, **self.pi_tensor_data})
            q_m, q_v, *_ = self.encoder(x, *args, **kwargs)
            output = q_m, q_v
        else:
            self.encoder.train(original_mode)
            raise NotImplementedError()
        self.encoder.train(original_mode)
        return output  # (K x L), (K x L)

    def log_prob(self, z: torch.Tensor) -> torch.Tensor:
        """Compute log probability of latent variables under the prior.

        Parameters
        ----------
        z
            Latent variables with shape (N, L).

        Returns
        -------
        torch.Tensor
            Log probability with shape (N, L).

        Notes
        -----
        The log probability is computed as:
        log p(z) = log Σ_k w_k * q(z|u_k)
        where the sum is computed using logsumexp for numerical stability.
        """
        # Mixture of gaussian computed on K x N x L
        z = z.unsqueeze(0)  # 1 x N x L

        # u->encoder->mean, var
        m_p, v_p = self.get_params()  # (K x L), (K x L)
        m_p = m_p.unsqueeze(1)  # K x 1 x L
        v_p = v_p.unsqueeze(1)  # K x 1 x L

        # mixing probabilities
        w = torch.nn.functional.softmax(self.w, dim=0)  # K x 1 x 1

        # sum of log_p across components weighted by w
        log_prob = Normal(m_p, v_p.sqrt()).log_prob(z) + torch.log(w)  # K x N x L
        log_prob = torch.logsumexp(log_prob, dim=0, keepdim=False)  # N x L

        return log_prob  # N x L

    def kl(self, qz: Normal) -> torch.Tensor:
        """Compute KL divergence using Monte Carlo estimation.

        Parameters
        ----------
        qz
            The posterior distribution.

        Returns
        -------
        torch.Tensor
            KL divergence between posterior and VampPrior.

        Notes
        -----
        The KL divergence is estimated using Monte Carlo sampling:
        KL(q||p) = E_q[log q(z) - log p(z)]
        where z ~ q(z) is sampled from the posterior.
        """
        assert isinstance(qz, Normal)
        z = qz.rsample()
        return qz.log_prob(z) - self.log_prob(z)

    def get_extra_state(self) -> dict[str, Any]:
        """Get extra state for serialization.

        Returns
        -------
        dict
            Extra state information.
        """
        return {
            "pi_aux_data": self.pi_aux_data,
            "input_type": self.input_type,
        }

    def set_extra_state(self, state: dict[str, Any]) -> None:
        """Set extra state from serialization.

        Parameters
        ----------
        state
            Extra state information.
        """
        self.pi_aux_data = state["pi_aux_data"]
        self.input_type = state["input_type"]


class GaussianMixtureModelPrior(Prior):
    """Gaussian Mixture Model prior distribution.

    This prior uses a mixture of Gaussian distributions with learnable
    parameters. Unlike VampPrior, the mixture components are not tied
    to the encoder network.

    Parameters
    ----------
    n_components
        Number of mixture components.
    n_latent
        Dimensionality of the latent space.
    data
        Initial means and variances as (means, variances).
    trainable_priors
        Whether the prior parameters should be trainable.

    Notes
    -----
    This prior learns the means, variances, and mixing weights of a
    Gaussian mixture model directly. It's simpler than VampPrior but
    may not capture data-specific structure as well.

    The prior is computed as:
    p(z) = Σ_k w_k * N(z|μ_k, σ_k²)
    where μ_k, σ_k², and w_k are learnable parameters.

    Examples
    --------
    >>> import torch
    >>> # Create GMM prior
    >>> prior = GaussianMixtureModelPrior(n_components=3, n_latent=5)
    >>> # Create a posterior distribution
    >>> qz = torch.distributions.Normal(torch.randn(10, 5), torch.ones(10, 5))
    >>> # Compute KL divergence
    >>> kl = prior.kl(qz)
    >>> print(kl.shape)  # torch.Size([10, 5])
    """

    # Based on VampPrior class

    def __init__(
        self,
        n_components: int,
        n_latent: int,
        data: tuple[torch.Tensor, torch.Tensor] | None = None,
        trainable_priors: bool = True,
    ) -> None:
        # Do we need python 2 compatibility?
        super().__init__()

        if data is None:
            p_m = torch.rand(n_components, n_latent)  # K x L
            p_v = torch.ones(n_components, n_latent)  # K x L
        else:
            p_m = data[0]
            p_v = data[1]
        self.p_m = torch.nn.Parameter(p_m, requires_grad=trainable_priors)
        self.p_v = torch.nn.Parameter(p_v, requires_grad=trainable_priors)

        # mixing weights
        self.w = torch.nn.Parameter(torch.zeros(self.p_m.shape[0], 1, 1))  # K x 1 x 1

    def log_prob(self, z: torch.Tensor) -> torch.Tensor:
        """Compute log probability of latent variables under the prior.

        Parameters
        ----------
        z
            Latent variables with shape (N, L).

        Returns
        -------
        torch.Tensor
            Log probability with shape (N, L).

        Notes
        -----
        The log probability is computed as:
        log p(z) = log Σ_k w_k * N(z|μ_k, σ_k²)
        where the sum is computed using logsumexp for numerical stability.
        """
        # Mixture of gaussian computed on K x N x L
        z = z.unsqueeze(0)  # 1 x N x L

        m_p = self.p_m.unsqueeze(1)  # K x 1 x L
        v_p = self.p_v.unsqueeze(1)  # K x 1 x L

        # mixing probabilities
        w = torch.nn.functional.softmax(self.w, dim=0)  # K x 1 x 1

        # sum of log_p across components weighted by w
        log_prob = Normal(m_p, v_p.sqrt()).log_prob(z) + torch.log(w)  # K x N x L
        log_prob = torch.logsumexp(log_prob, dim=0, keepdim=False)  # N x L

        return log_prob  # N x L

    def kl(self, qz: Normal) -> torch.Tensor:
        """Compute KL divergence using Monte Carlo estimation.

        Parameters
        ----------
        qz
            The posterior distribution.

        Returns
        -------
        torch.Tensor
            KL divergence between posterior and GMM prior.

        Notes
        -----
        The KL divergence is estimated using Monte Carlo sampling:
        KL(q||p) = E_q[log q(z) - log p(z)]
        where z ~ q(z) is sampled from the posterior.
        """
        assert isinstance(qz, Normal)
        z = qz.rsample()
        return qz.log_prob(z) - self.log_prob(z)
