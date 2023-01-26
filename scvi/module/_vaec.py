import numpy as np
import torch
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder, FCLayers

torch.backends.cudnn.benchmark = True


# Conditional VAE model
class VAEC(BaseModuleClass):
    """
    Conditional Variational auto-encoder model.

    This is an implementation of the CondSCVI model

    Parameters
    ----------
    n_input
        Number of input genes
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for the encoder and decoder neural network
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    """

    def __init__(
        self,
        n_input: int,
        n_labels: int = 0,
        n_hidden: Tunable[int] = 128,
        n_latent: Tunable[int] = 5,
        n_layers: Tunable[int] = 2,
        log_variational: bool = True,
        ct_weight: np.ndarray = None,
        dropout_rate: Tunable[float] = 0.05,
        **module_kwargs,
    ):
        super().__init__()
        self.dispersion = "gene"
        self.n_latent = n_latent
        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.dropout_rate = dropout_rate
        self.log_variational = log_variational
        self.gene_likelihood = "nb"
        self.latent_distribution = "normal"
        # Automatically deactivate if useless
        self.n_batch = 0
        self.n_labels = n_labels

        # gene dispersion
        self.px_r = torch.nn.Parameter(torch.randn(n_input))

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_cat_list=[n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            return_dist=True,
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=[n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
        )
        self.px_decoder = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_input), torch.nn.Softplus()
        )

        if ct_weight is not None:
            ct_weight = torch.tensor(ct_weight, dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)

    def _get_inference_input(self, tensors):
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]

        input_dict = dict(
            x=x,
            y=y,
        )
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]

        input_dict = {
            "z": z,
            "library": library,
            "y": y,
        }
        return input_dict

    @auto_move_data
    def inference(self, x, y, n_samples=1):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        library = x.sum(1).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)

        qz, z = self.z_encoder(x_, y)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand(
                (n_samples, library.size(0), library.size(1))
            )

        outputs = dict(z=z, qz=qz, library=library)
        return outputs

    @auto_move_data
    def generative(self, z, library, y):
        """Runs the generative model."""
        h = self.decoder(z, y)
        px_scale = self.px_decoder(h)
        px_rate = library * px_scale
        px = NegativeBinomial(px_rate, logits=self.px_r)
        return dict(px=px)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Loss computation."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        qz = inference_outputs["qz"]
        px = generative_outputs["px"]

        mean = torch.zeros_like(qz.loc)
        scale = torch.ones_like(qz.scale)

        kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)

        reconst_loss = -px.log_prob(x).sum(-1)
        scaling_factor = self.ct_weight[y.long()[:, 0]]
        loss = torch.mean(scaling_factor * (reconst_loss + kl_weight * kl_divergence_z))

        return LossOutput(
            loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence_z
        )

    @torch.inference_mode()
    def sample(
        self,
        tensors,
        n_samples=1,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Tensors dict
        n_samples
            Number of required samples for each cell

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        inference_kwargs = dict(n_samples=n_samples)
        generative_outputs = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )[1]

        px_r = generative_outputs["px_r"]
        px_rate = generative_outputs["px_rate"]

        dist = NegativeBinomial(px_rate, logits=px_r)
        if n_samples > 1:
            exprs = dist.sample().permute(
                [1, 2, 0]
            )  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample()

        return exprs.cpu()
