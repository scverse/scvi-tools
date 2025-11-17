import numpy as np
import torch
from torch.distributions import Categorical, Distribution, Independent, MixtureSameFamily, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.distributions import NegativeBinomial
from scvi.module.base import (
    BaseMinifiedModeModuleClass,
    EmbeddingModuleMixin,
    LossOutput,
    auto_move_data,
)
from scvi.nn import Encoder, FCLayers

from ._classifier import Classifier

torch.backends.cudnn.benchmark = True


class VAEC(EmbeddingModuleMixin, BaseMinifiedModeModuleClass):
    """Conditional Variational auto-encoder model.

    This is an implementation of the CondSCVI model

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches. If ``0``, no batch correction is performed.
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
    encode_covariates
        If ``True``, covariates are concatenated to gene expression prior to passing through
        the encoder(s). Else, only gene expression is used.
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
        n_fine_labels: int | None = None,
        n_hidden: int = 128,
        n_latent: int = 5,
        n_layers: int = 2,
        log_variational: bool = True,
        ct_weight: np.ndarray | None = None,
        dropout_rate: float = 0.05,
        encode_covariates: bool = False,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        linear_classifier: bool = True,
        prior: str = "normal",
        num_classes_mog: int | None = 10,
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
        self.n_fine_labels = n_fine_labels
        self.prior = prior
        self.num_classes_mog = num_classes_mog
        self.init_embedding(REGISTRY_KEYS.BATCH_KEY, n_batch, **{})

        if self.encode_covariates and self.n_batch < 1:
            raise ValueError("`n_batch` must be greater than 0 if `encode_covariates` is `True`.")

        batch_dim = self.get_embedding(REGISTRY_KEYS.BATCH_KEY).embedding_dim

        cat_list = [n_labels]
        n_input_encoder = n_input + batch_dim * encode_covariates

        # gene dispersion
        self.px_r = torch.nn.Parameter(torch.randn(n_input))

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        _extra_encoder_kwargs = {}
        self.z_encoder = Encoder(
            n_input_encoder,
            n_latent,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            return_dist=True,
            **(extra_encoder_kwargs or {}),
        )
        if n_fine_labels is not None:
            cls_parameters = {
                "n_layers": 0,
                "n_hidden": 0,
                "dropout_rate": dropout_rate,
                "logits": True,
            }
            # linear mapping from latent space to a coarse-celltype aware space
            self.linear_mapping = FCLayers(
                n_in=n_latent,
                n_out=n_hidden,
                n_cat_list=[n_labels],
                use_layer_norm=True,
                dropout_rate=0.0,
            )

            self.classifier = Classifier(
                n_hidden,
                n_labels=n_fine_labels,
                use_batch_norm=False,
                use_layer_norm=True,
                **cls_parameters,
            )
        else:
            self.classifier = None

        # decoder goes from n_latent-dimensional space to n_input-d data
        _extra_decoder_kwargs = {}
        n_input_decoder = n_latent + batch_dim
        self.decoder = FCLayers(
            n_in=n_input_decoder,
            n_out=n_hidden,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=True,
            use_batch_norm=False,
            use_layer_norm=True,
            **(extra_decoder_kwargs or {}),
        )
        self.px_decoder = torch.nn.Linear(n_hidden, n_input)
        self.per_ct_bias = torch.nn.Parameter(torch.zeros(n_labels, n_input))

        if ct_weight is not None:
            ct_weight = torch.tensor(ct_weight, dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)
        if self.prior == "mog":
            self.prior_means = torch.nn.Parameter(
                torch.randn([n_labels, self.num_classes_mog, n_latent])
            )
            self.prior_log_std = torch.nn.Parameter(
                torch.zeros([n_labels, self.num_classes_mog, n_latent]) - 2.0
            )
            self.prior_logits = torch.nn.Parameter(torch.zeros([n_labels, self.num_classes_mog]))

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
                raise NotImplementedError(f"Unknown minified-data type: {self.minified_data_type}")

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
            x_ = torch.log1p(x_)
        if self.encode_covariates:
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            encoder_input = torch.cat([x_, batch_rep], dim=-1)
        else:
            encoder_input = x_
        qz, z = self.z_encoder(encoder_input, y)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand((n_samples, library.size(0), library.size(1)))

        outputs = {"z": z, "qz": qz, "library": library}
        return outputs

    @auto_move_data
    def _cached_inference(
        self, qzm, qzv, observed_lib_size, y=None, batch_index=None, n_samples=1
    ):
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
            raise NotImplementedError(f"Unknown minified-data type: {self.minified_data_type}")
        outputs = {"z": z, "qz": qz, "library": library}
        return outputs

    @auto_move_data
    def classify(
        self,
        z: torch.Tensor,
        label_index: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward pass through the encoder and classifier.

        Parameters
        ----------
        z
            Tensor of shape ``(n_obs, n_latent)``.
        label_index
            Tensor of shape ``(n_obs,)`` denoting label indices.

        Returns
        -------
        Tensor of shape ``(n_obs, n_labels)`` denoting logit scores per label.
        """
        if len(label_index.shape) == 1:
            label_index = label_index.unsqueeze(1)
        classifier_latent = self.linear_mapping(z, label_index)
        w_y = self.classifier(classifier_latent)
        return w_y

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        y: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
    ) -> dict[str, Distribution]:
        """Runs the generative model."""
        batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
        # Handle case when z has an extra sample dimension (n_samples > 1)
        if z.ndim > batch_rep.ndim:
            batch_rep = batch_rep.unsqueeze(0).expand(z.shape[0], -1, -1)
        decoder_input = torch.cat([z, batch_rep], dim=-1)
        h = self.decoder(decoder_input, y, batch_index)
        px_scale = torch.nn.Softmax(dim=-1)(self.px_decoder(h) + self.per_ct_bias[y.ravel()])
        px_rate = library * px_scale
        px_r = torch.exp(self.px_r)
        px = NegativeBinomial(mu=px_rate, theta=px_r, scale=px_scale)
        return {"px": px}

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution],
        generative_outputs: dict[str, Distribution],
        kl_weight: float = 1.0,
        labelled_tensors: dict[str, torch.Tensor] | None = None,
        classification_ratio=5.0,
    ):
        """Loss computation."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY].ravel().long()
        qz = inference_outputs["qz"]
        px = generative_outputs["px"]
        fine_labels = (
            tensors["fine_labels"].ravel().long() if "fine_labels" in tensors.keys() else None
        )

        if self.prior == "mog":
            indexed_means = self.prior_means[y]
            indexed_log_std = self.prior_log_std[y]
            indexed_logits = self.prior_logits[y]
            cats = Categorical(logits=indexed_logits)
            normal_dists = Independent(
                Normal(indexed_means, torch.exp(indexed_log_std) + 1e-4),
                reinterpreted_batch_ndims=1,
            )
            prior = MixtureSameFamily(cats, normal_dists)
            u = qz.rsample(sample_shape=(30,))
            # (sample, n_obs, n_latent) -> (sample, n_obs,)
            kl_divergence_z = (qz.log_prob(u).sum(-1) - prior.log_prob(u)).mean(0)
        else:
            mean = torch.zeros_like(qz.loc)
            scale = torch.ones_like(qz.scale)
            kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)

        reconst_loss = -px.log_prob(x).sum(-1)
        scaling_factor = self.ct_weight[y]

        if self.classifier is not None:
            fine_labels = fine_labels.view(-1)
            logits = self.classify(
                qz.loc, label_index=tensors[REGISTRY_KEYS.LABELS_KEY]
            )  # (n_obs, n_labels)
            classification_loss_ = F.cross_entropy(logits, fine_labels, reduction="none")
            mask = fine_labels != self.n_fine_labels
            classification_loss = classification_ratio * torch.masked_select(
                classification_loss_, mask
            ).mean(0)

            loss = torch.mean(
                scaling_factor * (reconst_loss + classification_loss + kl_weight * kl_divergence_z)
            )

            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_divergence_z,
                classification_loss=classification_loss,
                logits=logits,
                true_labels=fine_labels,
            )
        else:
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
            exprs = dist.sample().permute([1, 2, 0])  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample()

        return exprs.cpu()
