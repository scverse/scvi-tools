import random
from typing import Any

import mlx.core as mx
import mlx.nn as nn

from scvi import REGISTRY_KEYS
from scvi.module.base import LossOutput


class MlxDense(nn.Module):
    """MLX dense layer with PyTorch-compatible initialization."""

    def __init__(self, in_features: int, out_features: int):
        super().__init__()
        scale = 1 / 3
        limit = scale * (6 / (in_features + out_features)) ** 0.5
        self.weight = mx.random.uniform(low=-limit, high=limit, shape=(in_features, out_features))
        self.bias = mx.zeros(out_features)

    def __call__(self, x):
        return x @ self.weight + self.bias


class MlxEncoder(nn.Module):
    """Encoder module for MLX VAE."""

    def __init__(self, n_input: int, n_latent: int, n_hidden: int, dropout_rate: float):
        super().__init__()
        self.dense1 = MlxDense(n_input, n_hidden)
        self.bn1 = nn.BatchNorm(n_hidden, momentum=0.9)
        self.dropout1 = nn.Dropout(dropout_rate)
        self.dense2 = MlxDense(n_hidden, n_hidden)
        self.bn2 = nn.BatchNorm(n_hidden, momentum=0.9)
        self.dropout2 = nn.Dropout(dropout_rate)
        self.mean_layer = MlxDense(n_hidden, n_latent)
        self.var_layer = MlxDense(n_hidden, n_latent)

    def __call__(self, x):
        x = mx.log1p(x)
        h = self.dense1(x)
        h = self.bn1(h)
        h = nn.relu(h)
        h = self.dropout1(h)
        h = self.dense2(h)
        h = self.bn2(h)
        h = nn.relu(h)
        h = self.dropout2(h)
        mean = self.mean_layer(h)
        log_var = self.var_layer(h)
        var = mx.exp(log_var)
        return mean, var


class MlxDecoder(nn.Module):
    """Decoder module for MLX VAE."""

    def __init__(
        self, n_input: int, n_hidden: int, n_batch: int, n_latent: int, dropout_rate: float = 0.0
    ):
        super().__init__()
        self.dense1 = MlxDense(n_latent, n_hidden)
        self.dense2 = MlxDense(n_batch, n_hidden)
        self.bn1 = nn.BatchNorm(n_hidden, momentum=0.9)
        self.dropout1 = nn.Dropout(dropout_rate)
        self.dense3 = MlxDense(n_hidden, n_hidden)
        self.dense4 = MlxDense(n_batch, n_hidden)
        self.bn2 = nn.BatchNorm(n_hidden, momentum=0.9)
        self.dropout2 = nn.Dropout(dropout_rate)
        self.dense5 = MlxDense(n_hidden, n_input)
        key = mx.random.key(random.randint(0, 2**31 - 1))
        self.disp = mx.random.normal(key=key, shape=(n_input,))

    def __call__(self, z, batch):
        h = self.dense1(z) + self.dense2(batch)
        h = self.bn1(h)
        h = nn.relu(h)
        h = self.dropout1(h)
        h = self.dense3(h) + self.dense4(batch)
        h = self.bn2(h)
        h = nn.relu(h)
        h = self.dropout2(h)
        rho_unnorm = self.dense5(h)
        return rho_unnorm, self.disp


# TODO: WE CAN INHERIT BaseModuleClass instead of nn.Module AND MODEL SAVE/LOAD WILL BE OK BUT
# WITH SOME ISSUES WITH GETTING THE TRAINING LOSS (otherwise use BaseModuleClass)
class MlxVAE(nn.Module):
    """Variational Autoencoder model using the MLX framework."""

    def __init__(
        self,
        n_input: int,
        n_batch: int,
        n_hidden: int = 128,
        n_latent: int = 30,
        dropout_rate: float = 0.0,
        n_layers: int = 1,
        gene_likelihood: str = "nb",
        eps: float = 1e-8,
    ):
        super().__init__()
        self.n_input = n_input
        self.n_batch = n_batch
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.dropout_rate = dropout_rate
        self.n_layers = n_layers
        self.gene_likelihood = gene_likelihood
        self.eps = eps
        self.encoder = MlxEncoder(n_input, n_latent, n_hidden, dropout_rate)
        self.decoder = MlxDecoder(n_input, n_hidden, n_batch, n_latent, dropout_rate)

    def state_dict(self) -> dict[str, Any]:
        """Return the model parameters as a state dictionary for saving.

        This method provides compatibility with the scvi-tools save/load mechanism
        by converting MLX parameters to a format that can be serialized.
        """
        import mlx.core as mx

        # Get all parameters as a nested dictionary and convert to numpy arrays
        params = self.parameters()

        def convert_to_numpy(obj):
            if isinstance(obj, mx.array):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_to_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_to_numpy(item) for item in obj]
            else:
                return obj

        return convert_to_numpy(params)

    def load_state_dict(self, state_dict: dict[str, Any]) -> None:
        """Load parameters from a state dictionary.

        Parameters
        ----------
        state_dict
            Dictionary containing model parameters.
        """
        import mlx.core as mx

        def convert_to_mlx(obj):
            if isinstance(obj, list):
                # Check if it's a list of numbers (i.e., an array)
                if obj and isinstance(obj[0], (int, float)):
                    return mx.array(obj)
                elif obj and isinstance(obj[0], list):
                    return mx.array(obj)
                else:
                    return [convert_to_mlx(item) for item in obj]
            elif isinstance(obj, dict):
                return {k: convert_to_mlx(v) for k, v in obj.items()}
            else:
                return obj

        mlx_params = convert_to_mlx(state_dict)
        self.update(mlx_params)

    def on_load(self, model, **kwargs) -> None:
        """Callback function run after loading a saved model.

        This method is called during model loading to perform any necessary
        post-load initialization. For MLX models, this is a no-op since
        MLX doesn't use Pyro or require special parameter store handling.

        Parameters
        ----------
        model
            The model instance being loaded.
        **kwargs
            Additional keyword arguments (e.g., pyro_param_store for Pyro models).
        """
        pass

    def _get_inference_input(self, tensors: dict[str, mx.array]) -> dict[str, mx.array]:
        """Get input for the inference model."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        return {"x": x}

    def _get_generative_input(
        self,
        tensors: dict[str, mx.array],
        inference_outputs: dict[str, mx.array],
    ) -> dict[str, mx.array]:
        """Get input for the generative model."""
        z = inference_outputs["z"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        return {
            "z": z,
            "batch_index": batch_index,
        }

    def train(self, mode=True):
        """Recursively set child modules to training mode."""
        super().train(mode)
        for module in self.children():
            if isinstance(module, nn.Module):
                module.train(mode)
        return self

    def eval(self):
        """Recursively set child modules to evaluation mode."""
        super().eval()
        for module in self.children():
            if isinstance(module, nn.Module):
                module.eval()
        return self

    def inference(self, x: mx.array, n_samples: int = 1) -> dict[str, Any]:
        mean, var = self.encoder(x)
        stddev = mx.sqrt(mx.clip(var, 1e-8, 50.0)) + self.eps
        key = mx.random.key(random.randint(0, 2**31 - 1))
        shape = (n_samples,) + mean.shape
        eps = mx.random.normal(key=key, shape=shape)
        z = mean + stddev * eps
        mx.eval(z)
        return {"mean": mean, "var": var, "z": z}

    def generative(self, z, batch_index) -> dict[str, Any]:
        batch = mx.zeros((batch_index.shape[0], self.n_batch))
        rows = mx.arange(batch_index.reshape(-1).shape[0])
        batch = batch.at[rows, batch_index.reshape(-1)].add(1.0)
        rho_unnorm, disp = self.decoder(z, batch)
        rho_unnorm = mx.clip(rho_unnorm, -50, 50)
        rho_exp = mx.exp(rho_unnorm)
        rho = rho_exp / mx.sum(rho_exp, axis=-1, keepdims=True)
        return {"rho": rho, "disp": disp}

    def loss(
        self, tensors, inference_outputs, generative_outputs, kl_weight: float = 1.0
    ) -> LossOutput:
        x = tensors[REGISTRY_KEYS.X_KEY]
        disp = generative_outputs["disp"]
        mean = inference_outputs["mean"]
        var = inference_outputs["var"]
        rho = generative_outputs["rho"]
        total_count = mx.sum(x, axis=-1, keepdims=True)
        mu = total_count * rho

        eps = 1e-10
        # Clip dispersion to prevent numerical instability
        disp = mx.clip(disp, eps, 50.0)
        # Clip mu to prevent extreme values
        mu = mx.clip(mu, eps, 1e6)

        log_theta_mu_eps = mx.log(disp + mu + eps)
        log_theta_eps = mx.log(disp + eps)
        log_mu_eps = mx.log(mu + eps)
        log_prob = x * (log_mu_eps - log_theta_mu_eps) + disp * (log_theta_eps - log_theta_mu_eps)

        # Clip arguments to log1p to prevent overflow
        log_theta_mu_eps_clipped = mx.clip(log_theta_mu_eps, -50, 50)
        log_theta_eps_clipped = mx.clip(log_theta_eps, -50, 50)
        log_prob += mx.log1p(mx.exp(log_theta_mu_eps_clipped)) - mx.log1p(
            mx.exp(log_theta_eps_clipped)
        )

        reconst_loss = -mx.sum(log_prob, axis=-1)
        var = mx.clip(var, eps, 50.0)
        kl_divergence = 0.5 * mx.sum(
            var + mx.square(mean) - 1.0 - mx.log(mx.maximum(var, eps) + eps), axis=-1
        )
        weighted_kl = kl_weight * kl_divergence
        loss = mx.mean(reconst_loss + weighted_kl)
        return LossOutput(loss=-loss, reconstruction_loss=-reconst_loss, kl_local=-kl_divergence)

    def __call__(
        self, tensors: dict[str, mx.array], kl_weight: float = 1.0
    ) -> tuple[dict, dict, LossOutput]:
        inference_inputs = self._get_inference_input(tensors)
        inference_outputs = self.inference(**inference_inputs)
        generative_inputs = self._get_generative_input(tensors, inference_outputs)
        generative_outputs = self.generative(**generative_inputs)
        loss_output = self.loss(tensors, inference_outputs, generative_outputs, kl_weight)
        return inference_outputs, generative_outputs, loss_output
