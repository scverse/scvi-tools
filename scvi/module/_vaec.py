from collections.abc import Iterable
from typing import Optional

import numpy as np
import torch
from torch.distributions import Categorical, Independent, MixtureSameFamily, Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi._types import Tunable
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.distributions import NegativeBinomial
from scvi.module.base import BaseMinifiedModeModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder, FCLayers

torch.backends.cudnn.benchmark = True


# Conditional VAE model
class VAEC(BaseMinifiedModeModuleClass):
    """Conditional Variational auto-encoder model.

    This is an implementation of the CondSCVI model

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    ct_weight
        Multiplicative weight for cell type specific latent space.
    dropout_rate
        Dropout rate for the encoder and decoder neural network.
    extra_encoder_kwargs
        Keyword arguments passed into :class:`~scvi.nn.Encoder`.
    extra_decoder_kwargs
        Keyword arguments passed into :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: Tunable[int] = 128,
        n_latent: Tunable[int] = 5,
        n_layers: Tunable[int] = 2,
        log_variational: bool = True,
        ct_weight: np.ndarray = None,
        dropout_rate: Tunable[float] = 0.05,
        encode_covariates: bool = False,
        extra_encoder_kwargs: Optional[dict] = None,
        extra_decoder_kwargs: Optional[dict] = None,
        prior: str = 'normal',
        df_ct_id_dict: dict = None,
        num_classes_mog: Optional[int] = 10,
    ):
        super().__init__()
        self.dispersion = "gene"
        self.n_latent = n_latent
        self.n_layers = n_layers
        self.n_hidden = n_hidden
        self.dropout_rate = dropout_rate
        self.encode_covariates = encode_covariates
        self.log_variational = log_variational
        self.gene_likelihood = "nb"
        self.latent_distribution = "normal"
        # Automatically deactivate if useless
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.prior = prior
        if df_ct_id_dict is not None:
            self.num_classes_mog = max([v[2] for v in df_ct_id_dict.values()]) + 1
            mapping_mog = torch.tensor([v[2] for _, v in sorted(df_ct_id_dict.items())])
            self.register_buffer("mapping_mog", mapping_mog)
        else:
            self.num_classes_mog = num_classes_mog
        cat_list = [n_labels, n_batch]
        encoder_cat_list = cat_list if self.encode_covariates else [n_labels]

        # gene dispersion
        self.px_r = torch.nn.Parameter(torch.randn(n_input))

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        _extra_encoder_kwargs = {}
        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_cat_list=encoder_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            return_dist=True,
            **_extra_encoder_kwargs,
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        _extra_decoder_kwargs = {}
        self.decoder = FCLayers(
            n_in=n_latent,
            n_out=n_hidden,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            **_extra_decoder_kwargs,
        )
        self.px_decoder = torch.nn.Sequential(
            torch.nn.Linear(n_hidden, n_input), torch.nn.Softplus()
        )

        if ct_weight is not None:
            ct_weight = torch.tensor(ct_weight, dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)
        if self.prior=='mog':
            self.prior_means = torch.nn.Parameter(
                    0.01 * torch.randn([n_labels, self.num_classes_mog, n_latent]))
            self.prior_log_scales = torch.nn.Parameter(
                    torch.zeros([n_labels, self.num_classes_mog, n_latent]))
            self.prior_logits = torch.nn.Parameter(
                    torch.zeros([n_labels, self.num_classes_mog]))

    def _get_inference_input(self, tensors):
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        if self.minified_data_type is None:
            x = tensors[REGISTRY_KEYS.X_KEY]
            input_dict = {
                "x": x,
                "y": y,
                "batch_index": batch_index,
            }
        else:
            if ADATA_MINIFY_TYPE.__contains__(self.minified_data_type):
                qzm = tensors[REGISTRY_KEYS.LATENT_QZM_KEY]
                qzv = tensors[REGISTRY_KEYS.LATENT_QZV_KEY]
                observed_lib_size = tensors[REGISTRY_KEYS.OBSERVED_LIB_SIZE]
                input_dict = {
                    "qzm": qzm,
                    "qzv": qzv,
                    "observed_lib_size": observed_lib_size,
                    "y": y,
                    "batch_index": batch_index,
                }
            else:
                raise NotImplementedError(
                    f"Unknown minified-data type: {self.minified_data_type}"
                )

        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            "z": z,
            "library": library,
            "y": y,
            "batch_index": batch_index,
        }
        return input_dict

    @auto_move_data
    def _regular_inference(self, x, y, batch_index, n_samples=1):
        """High level inference method.

        Runs the inference (encoder) model.
        """
        x_ = x
        library = x.sum(1).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)
        if self.encode_covariates:
            categorical_input = [y, batch_index]
        else:
            categorical_input = [y]
        qz, z = self.z_encoder(x_, *categorical_input)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand(
                (n_samples, library.size(0), library.size(1))
            )

        outputs = {"z": z, "qz": qz, "library": library}
        return outputs

    @auto_move_data
    def _cached_inference(self, qzm, qzv, observed_lib_size, n_samples=1):
        if ADATA_MINIFY_TYPE.__contains__(self.minified_data_type):
            qz = Normal(qzm, qzv.sqrt())
            # use dist.sample() rather than rsample because we aren't optimizing the z here
            untran_z = qz.sample() if n_samples == 1 else qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = observed_lib_size
            if n_samples > 1:
                library = library.unsqueeze(0).expand(
                    (n_samples, library.size(0), library.size(1))
                )
        else:
            raise NotImplementedError(
                f"Unknown minified-data type: {self.minified_data_type}"
            )
        outputs = {"z": z, "qz": qz, "library": library}
        return outputs

    @auto_move_data
    def generative(self, z, library, y, batch_index):
        """Runs the generative model."""
        h = self.decoder(z, y, batch_index)
        px_scale = self.px_decoder(h)
        px_rate = library * px_scale
        px = NegativeBinomial(px_rate, logits=self.px_r)
        return {"px": px, "px_scale": px_scale}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Loss computation."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY].ravel().long()
        qz = inference_outputs["qz"]
        px = generative_outputs["px"]
        fine_celltypes = tensors['fine_labels'].ravel().long() if 'fine_labels' in tensors.keys() else None

        if self.prior == "mog":
            indexed_means = self.prior_means[y]
            indexed_log_scales = self.prior_log_scales[y]
            indexed_logits = self.prior_logits[y]

            # Assigns zero meaning equal weight to all unlabeled cells. Otherwise biases to sample from respective MoG.
            if fine_celltypes is not None:
                logits_input = torch.nn.functional.one_hot(
                    self.mapping_mog[fine_celltypes], self.num_classes_mog)
                cats = Categorical(logits=10*logits_input + indexed_logits)
            else:
                cats = Categorical(logits=indexed_logits)
            normal_dists = torch.distributions.Independent(
                Normal(
                    indexed_means,
                    torch.exp(indexed_log_scales) + 1e-4
                ),
                reinterpreted_batch_ndims=1
            )
            prior = MixtureSameFamily(cats, normal_dists)
            u = qz.rsample(sample_shape=(30,))
            # (sample, n_obs, n_latent) -> (sample, n_obs,)
            kl_divergence_z = - (prior.log_prob(u) - qz.log_prob(u).sum(-1)).mean(0)
        else:
            mean = torch.zeros_like(qz.loc)
            scale = torch.ones_like(qz.scale)
            kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)

        reconst_loss = -px.log_prob(x).sum(-1)
        scaling_factor = self.ct_weight[y]
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
        r"""Generate observation samples from the posterior predictive distribution.

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
        inference_kwargs = {"n_samples": n_samples}
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
