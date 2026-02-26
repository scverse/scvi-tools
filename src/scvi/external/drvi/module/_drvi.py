from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import torch
from torch.distributions import Normal, kl_divergence

from scvi import REGISTRY_KEYS
from scvi.external.drvi.module._constants import MODULE_KEYS
from scvi.external.drvi.nn import DecoderDRVI, Encoder
from scvi.external.drvi.nn_modules.embedding import MultiEmbedding
from scvi.external.drvi.nn_modules.noise_model import (
    LogNegativeBinomialNoiseModel,
    NegativeBinomialNoiseModel,
    NormalNoiseModel,
    PoissonNoiseModel,
)
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence
    from typing import Any, Literal

    from scvi.external.drvi.nn_modules.layer.factory import LayerFactory

    TensorDict = dict[str, torch.Tensor]


class DRVIModule(BaseModuleClass):
    """DRVI (Disentangled Representation Variational Inference) pytorch module.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_latent
        Dimensionality of the latent space.
    n_split_latent
        Number of splits in the latent space. -1 means split all dimensions
        (n_split_latent=n_latent).
    split_aggregation
        How to aggregate splits in the last layer of the decoder.
    split_method
        How to make splits:
        - "split" : Split the latent space
        - "power" : Transform the latent space to n_split vectors of size n_latent
        - "power@Z" : Transform the latent space to n_split vectors of size n_latent Z
        - "split_map" : Split the latent space then map each to latent space using
          unique transformations
        - "split_map@Z" : Split the latent space then map each to vector of size Z
          using unique transformations
        - "split_diag" : Simple diagonal splitting
    decoder_reuse_weights
        Where to reuse the weights of the decoder layers when using splitting.
        Possible values are 'everywhere', 'last', 'intermediate', 'nowhere',
        'not_first'. Defaults to "everywhere".
    encoder_dims
        Number of nodes in hidden layers of the encoder.
    decoder_dims
        Number of nodes in hidden layers of the decoder.
    n_cats_per_cov
        Number of categories for each categorical covariate.
    n_continuous_cov
        Number of continuous covariates.
    encode_covariates
        Whether to concatenate covariates to expression in encoder.
    deeply_inject_covariates
        Whether to concatenate covariates into output of hidden layers in encoder/decoder.
        This option only applies when `n_layers` >= 1. The covariates are concatenated
        to the input of subsequent hidden layers.
    categorical_covariate_dims
        Embedding dimension of covariate keys if applicable.
    covariate_modeling_strategy
        The strategy model takes to remove covariates.
    use_batch_norm
        Whether to use batch norm in layers.
    affine_batch_norm
        Whether to use affine batch norm in layers.
    use_layer_norm
        Whether to use layer norm in layers.
    fill_in_the_blanks_ratio
        Ratio for fill-in-the-blanks training.
    reconstruction_strategy
        Strategy for reconstruction.
        - "dense" : Reconstruct all features.
        - "random_batch@M" : Reconstruct M random features for each batch.
    input_dropout_rate
        Dropout rate to apply to the input.
    encoder_dropout_rate
        Dropout rate to apply to each of the encoder hidden layers.
    decoder_dropout_rate
        Dropout rate to apply to each of the decoder hidden layers.
    gene_likelihood
        Gene likelihood model. Options include:
        - "poisson" : Poisson distributions
        - "nb" : Negative binomial distributions
        - "pnb": Log negative binomial distributions
        - "normal" : Normal distributions
        - "normal_unit_var" : Normal distributions with unit variance
    prior
        Prior model.
    var_activation
        The activation function to ensure positivity of the variational distribution.
        Options include "exp", "pow2", "2sig" or a custom callable.
    mean_activation
        The activation function at the end of mean encoder.
        Options include "identity", "relu", "leaky_relu", "leaky_relu_{slope}",
        "elu", "elu_{min_value}" or a custom callable.
    encoder_layer_factory
        A layer Factory instance for building encoder layers.
    decoder_layer_factory
        A layer Factory instance for building decoder layers.
    extra_encoder_kwargs
        Extra keyword arguments passed into encoder.
    extra_decoder_kwargs
        Extra keyword arguments passed into decoder.
    """

    def __init__(
        self,
        n_input: int,
        n_latent: int = 32,
        n_split_latent: int | None = -1,
        split_aggregation: Literal["sum", "logsumexp", "max"] = "logsumexp",
        split_method: Literal["split", "power", "split_map", "split_diag"] | str = "split_map",
        decoder_reuse_weights: Literal[
            "everywhere", "last", "intermediate", "nowhere", "not_first"
        ] = "everywhere",
        encoder_dims: Sequence[int] = (128, 128),
        decoder_dims: Sequence[int] = (128, 128),
        n_cats_per_cov: Iterable[int] | None = (),
        n_continuous_cov: int = 0,
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = False,
        categorical_covariate_dims: Sequence[int] = (),
        covariate_modeling_strategy: Literal[
            "one_hot",
            "emb",
            "emb_shared",
            "one_hot_linear",
            "emb_linear",
            "emb_shared_linear",
        ] = "one_hot",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        affine_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        fill_in_the_blanks_ratio: float = 0.0,
        reconstruction_strategy: str = "dense",
        input_dropout_rate: float = 0.0,
        encoder_dropout_rate: float = 0.1,
        decoder_dropout_rate: float = 0.0,
        gene_likelihood: Literal["normal", "poisson", "nb", "pnb"] = "pnb",
        dispersion: Literal["gene", "gene-batch", "gene-cell"] = "gene",
        prior: Literal["normal"] = "normal",
        var_activation: Callable | Literal["exp", "pow2", "2sig"] = "exp",
        mean_activation: Callable | str = "identity",
        encoder_layer_factory: LayerFactory | None = None,
        decoder_layer_factory: LayerFactory | None = None,
        extra_encoder_kwargs: dict[str, Any] | None = None,
        extra_decoder_kwargs: dict[str, Any] | None = None,
    ) -> None:
        super().__init__()
        self.n_latent = n_latent
        self.encoder_dims = encoder_dims
        self.decoder_dims = decoder_dims
        if n_split_latent is None or n_split_latent == -1:
            n_split_latent = n_latent
        self.n_split_latent = n_split_latent
        self.split_aggregation = split_aggregation
        self.latent_distribution = "normal"
        self.gene_likelihood = gene_likelihood

        self.encode_covariates = encode_covariates
        self.deeply_inject_covariates = deeply_inject_covariates

        self.gene_likelihood_module = self._construct_gene_likelihood_module(
            gene_likelihood, dispersion
        )
        self.fill_in_the_blanks_ratio = fill_in_the_blanks_ratio
        self.reconstruction_strategy = reconstruction_strategy

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        affine_batch_norm_encoder = affine_batch_norm == "encoder" or affine_batch_norm == "both"
        affine_batch_norm_decoder = affine_batch_norm == "decoder" or affine_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        assert covariate_modeling_strategy in [
            "one_hot",
            "emb",
            "emb_shared",
            "one_hot_linear",
            "emb_linear",
            "emb_shared_linear",
        ]
        if (
            covariate_modeling_strategy in ["emb_shared", "emb_shared_linear"]
            and len(n_cats_per_cov) > 0
        ):
            self.shared_covariate_emb = MultiEmbedding(
                n_cats_per_cov, categorical_covariate_dims, init_method="normal", max_norm=1.0
            )
        else:
            self.register_module("shared_covariate_emb", None)

        self.z_encoder = Encoder(
            n_input,
            n_latent,
            layers_dim=encoder_dims,
            input_dropout_rate=input_dropout_rate,
            dropout_rate=encoder_dropout_rate,
            n_cat_list=n_cats_per_cov if self.encode_covariates else [],
            n_continuous_cov=n_continuous_cov if self.encode_covariates else 0,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_encoder,
            affine_batch_norm=affine_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            var_activation=var_activation,
            mean_activation=mean_activation,
            layer_factory=encoder_layer_factory,
            covariate_modeling_strategy=covariate_modeling_strategy,
            categorical_covariate_dims=categorical_covariate_dims
            if self.encode_covariates
            else [],
            **(extra_encoder_kwargs or {}),
        )

        self.decoder = DecoderDRVI(
            n_latent,
            n_input,
            n_split=n_split_latent,
            split_aggregation=split_aggregation,
            split_method=split_method,
            reuse_weights=decoder_reuse_weights,
            gene_likelihood_module=self.gene_likelihood_module,
            layers_dim=decoder_dims,
            dropout_rate=decoder_dropout_rate,
            n_cat_list=n_cats_per_cov,
            n_continuous_cov=n_continuous_cov,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            affine_batch_norm=affine_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            layer_factory=decoder_layer_factory,
            covariate_modeling_strategy=covariate_modeling_strategy,
            categorical_covariate_dims=categorical_covariate_dims,
            **(extra_decoder_kwargs or {}),
        )

        self.prior = prior
        self.inspect_mode = False

    @property
    def fully_deterministic(self) -> bool:
        return self.z_encoder.fully_deterministic

    @fully_deterministic.setter
    def fully_deterministic(self, value: bool) -> None:
        self.z_encoder.fully_deterministic = value

    def _construct_gene_likelihood_module(self, gene_likelihood: str, dispersion: str) -> Any:
        """Construct the gene likelihood module based on the specified type.

        Parameters
        ----------
        gene_likelihood
            Type of gene likelihood model to construct.
        dispersion
            Dispersion parameter modeling strategy. Only used for relevant likelihoods.

        Returns
        -------
        object
            Constructed gene likelihood module.

        Raises
        ------
        NotImplementedError
            If the gene likelihood type is not supported.
        """
        if gene_likelihood == "normal_unit_var":
            return NormalNoiseModel(model_var="fixed=1")
        elif gene_likelihood == "normal":
            return NormalNoiseModel(model_var=dispersion)
        elif gene_likelihood == "poisson":
            return PoissonNoiseModel(mean_transformation="softmax", library_normalization="none")
        elif gene_likelihood in ["nb"]:
            return NegativeBinomialNoiseModel(
                dispersion=dispersion, mean_transformation="softmax", library_normalization="none"
            )
        elif gene_likelihood in ["pnb"]:
            return LogNegativeBinomialNoiseModel(
                dispersion=dispersion, mean_transformation="softmax", library_normalization="none"
            )
        else:
            raise NotImplementedError()

    def _get_inference_input(self, tensors: TensorDict) -> dict[str, torch.Tensor | None]:
        """Parse the dictionary to get appropriate args.

        Parameters
        ----------
        tensors
            Dictionary containing tensor data.

        Returns
        -------
        dict
            Dictionary with parsed input data.
        """
        x = tensors[REGISTRY_KEYS.X_KEY]

        batch_index = tensors.get(REGISTRY_KEYS.BATCH_KEY)
        cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY)
        cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY)

        if batch_index is not None:
            if cat_covs is not None:
                cat_covs = torch.cat([batch_index, cat_covs], dim=1)
            else:
                cat_covs = batch_index

        input_dict = {
            MODULE_KEYS.X_KEY: x,
            MODULE_KEYS.CONT_COVS_KEY: cont_covs,
            MODULE_KEYS.CAT_COVS_KEY: cat_covs,
        }
        return input_dict

    def _input_pre_processing(
        self,
        x: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
    ) -> dict[str, Any]:
        """Pre-process input data for the model.

        Parameters
        ----------
        x
            Input tensor data.
        cont_covs
            Continuous covariates.
        cat_covs
            Categorical covariates.

        Returns
        -------
        dict
            Dictionary containing pre-processed input data.
        """
        # log the input to the variational distribution for numerical stability
        x_ = self.gene_likelihood_module.initial_transformation(x)

        encoder_input = x_

        return {
            MODULE_KEYS.X_KEY: encoder_input,
            MODULE_KEYS.CAT_COVS_KEY: cat_covs if self.encode_covariates else None,
            MODULE_KEYS.CONT_COVS_KEY: cont_covs if self.encode_covariates else None,
        }

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, Any]:
        """High level inference method.

        Runs the inference (encoder) model.

        Parameters
        ----------
        x
            Input tensor data.
        cont_covs
            Continuous covariates.
        cat_covs
            Categorical covariates.
        n_samples
            Number of samples to generate.

        Returns
        -------
        dict
            Dictionary containing inference outputs including latent variables.
        """
        pre_processed_input = self._input_pre_processing(x, cont_covs, cat_covs).copy()
        x_ = pre_processed_input[MODULE_KEYS.X_KEY]
        library = self._get_library_size(x_original=x, reconstruction_indices=None)

        # Mask if needed
        if self.fill_in_the_blanks_ratio > 0.0 and self.training:
            assert cont_covs is None  # We do not consider cont_cov here
            x_mask = torch.where(torch.rand_like(x_) >= self.fill_in_the_blanks_ratio, 1.0, 0.0)
            x_ = x_ * x_mask
        else:
            x_mask = None

        # Prepare shared emb
        if self.shared_covariate_emb is not None and self.encode_covariates:
            pre_processed_input[MODULE_KEYS.CAT_COVS_KEY] = self.shared_covariate_emb(
                pre_processed_input[MODULE_KEYS.CAT_COVS_KEY].int()
            )

        # get variational parameters via the encoder networks
        qz_m, qz_v, z = self.z_encoder(
            x_,
            cat_full_tensor=pre_processed_input[MODULE_KEYS.CAT_COVS_KEY],
            cont_full_tensor=pre_processed_input[MODULE_KEYS.CONT_COVS_KEY],
        )

        outputs = {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.QZM_KEY: qz_m,
            MODULE_KEYS.QZV_KEY: qz_v,
            MODULE_KEYS.QL_KEY: None,  # We do not model library size
            MODULE_KEYS.LIBRARY_KEY: library,
            MODULE_KEYS.X_MASK_KEY: x_mask,
            MODULE_KEYS.N_SAMPLES_KEY: n_samples,
        }

        if n_samples > 1:
            for key in [
                MODULE_KEYS.Z_KEY,
                MODULE_KEYS.QZM_KEY,
                MODULE_KEYS.QZV_KEY,
                MODULE_KEYS.LIBRARY_KEY,
                MODULE_KEYS.X_MASK_KEY,
            ]:
                if outputs[key] is None:
                    continue
                assert outputs[key].shape[0] == z.shape[0]
                outputs[key] = (
                    outputs[key].unsqueeze(0).repeat(n_samples, *([1] * outputs[key].ndim))
                )
            outputs[MODULE_KEYS.Z_KEY] = Normal(
                outputs[MODULE_KEYS.QZM_KEY], outputs[MODULE_KEYS.QZV_KEY].sqrt()
            ).rsample()

        return outputs

    def _get_reconstruction_indices(self, tensors: TensorDict) -> None | torch.Tensor:
        # We also reconstruct a fraction in validation set
        if not self.training:
            return None
        if self.reconstruction_strategy == "dense":
            return None
        elif self.reconstruction_strategy.startswith("random_batch@"):
            x = tensors[REGISTRY_KEYS.X_KEY]
            n_random_features = int(self.reconstruction_strategy.split("@")[1])
            random_indices = torch.randperm(x.shape[1])[:n_random_features]
            return random_indices
        else:
            raise NotImplementedError(
                f"Reconstruction strategy {self.reconstruction_strategy} not implemented."
            )

    def _get_library_size(
        self, x_original: TensorDict, reconstruction_indices: torch.Tensor | None = None
    ) -> torch.Tensor:
        # Note: this is different from scvi implementation of library size that is log transformed
        # All our noise models accept non-normalized library to work
        if reconstruction_indices is None:
            return x_original.sum(1)
        elif reconstruction_indices.dim() == 1:
            return x_original[:, reconstruction_indices.to(x_original.device)].sum(1)
        else:
            raise NotImplementedError(
                f"Reconstruction indices {reconstruction_indices} not implemented."
            )

    def _get_generative_input(
        self,
        tensors: TensorDict,
        inference_outputs: dict[str, Any],
        transform_batch: int | None = None,
        library_to_inject: torch.Tensor | None = None,
    ) -> dict[str, Any]:
        """Prepare input for the generative model.

        Parameters
        ----------
        tensors
            Dictionary containing tensor data.
        inference_outputs
            Outputs from the inference step.
        library_to_inject
            Library size to inject (it should not be log transformed. Just x.sum(1)).
        transform_batch
            Batch to condition on.

        Returns
        -------
        dict
            Dictionary containing input for generative model.
        """
        z = inference_outputs[MODULE_KEYS.Z_KEY]

        batch_index = tensors.get(REGISTRY_KEYS.BATCH_KEY)
        cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY)
        cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY)

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        if batch_index is not None:
            if cat_covs is not None:
                cat_covs = torch.cat([batch_index, cat_covs], dim=1)
            else:
                cat_covs = batch_index

        n_samples = inference_outputs.get(MODULE_KEYS.N_SAMPLES_KEY, 1)

        reconstruction_indices = self._get_reconstruction_indices(tensors)

        # Set library size
        if library_to_inject is not None:
            library = library_to_inject
            assert reconstruction_indices is None
        elif (
            reconstruction_indices is not None
        ):  # Override library size as we do not decode everything
            library = self._get_library_size(tensors[REGISTRY_KEYS.X_KEY], reconstruction_indices)
        elif MODULE_KEYS.LIBRARY_KEY in inference_outputs:
            library = inference_outputs[MODULE_KEYS.LIBRARY_KEY]
        else:
            raise ValueError("Library size not found in inference outputs.")

        input_dict = {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.LIBRARY_KEY: library,
            MODULE_KEYS.CONT_COVS_KEY: cont_covs,
            MODULE_KEYS.CAT_COVS_KEY: cat_covs,
            MODULE_KEYS.N_SAMPLES_KEY: n_samples,
            MODULE_KEYS.RECONSTRUCTION_INDICES: reconstruction_indices,
        }

        if n_samples > 1:
            # Repeat the covariates for each sample
            for key in [MODULE_KEYS.CAT_COVS_KEY, MODULE_KEYS.CONT_COVS_KEY]:
                if input_dict[key] is None:
                    continue
                input_dict[key] = input_dict[key].repeat(
                    n_samples, *([1] * (input_dict[key].ndim - 1))
                )
            # Combine samples and batch dimensions into the first dimension
            for key in [MODULE_KEYS.Z_KEY, MODULE_KEYS.LIBRARY_KEY]:
                if input_dict[key] is not None:
                    input_dict[key] = input_dict[key].flatten(0, 1)
        return input_dict

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        cat_covs: torch.Tensor | None = None,
        cont_covs: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
        reconstruction_indices: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, Any]:
        """Runs the generative model.

        Parameters
        ----------
        z
            Latent variables.
        library
            Library size information.
        cont_covs
            Continuous covariates.
        cat_covs
            Categorical covariates.
        transform_batch
            Batch to condition on. Currently not used but required for RNASeqMixin compatibility.
        reconstruction_indices
            Indices of features to reconstruct.

        Returns
        -------
        dict
            Dictionary containing generative model outputs.
        """
        # Parameter transform_batch is not used!
        # But, we keep it here since _rna_mixin.py checks module.generative
        # to include this as a parameter!

        if self.shared_covariate_emb is not None:
            cat_covs = self.shared_covariate_emb(cat_covs.int())
        # form the likelihood
        px, params, original_params = self.decoder(
            z,
            cat_full_tensor=cat_covs,
            cont_full_tensor=cont_covs,
            library=library,
            reconstruction_indices=reconstruction_indices,
            return_original_params=self.inspect_mode,
        )

        if n_samples > 1:
            n_batch = z.shape[0] // n_samples
            for key, value in params.items():
                value = value.reshape(
                    n_samples, n_batch, *value.shape[1:]
                )  # Shape : (n_batch, n_samples, n_genes)
                params[key] = value

            library = library.reshape(n_samples, n_batch, *library.shape[1:])
            px = self.gene_likelihood_module.dist(parameters=params, lib_y=library)

        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PX_PARAMS_KEY: params,
            MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY: original_params,
            MODULE_KEYS.RECONSTRUCTION_INDICES: reconstruction_indices,
        }

    def _get_reconstruction_loss(
        self,
        px: torch.distributions.Distribution,
        x: torch.Tensor,
        x_mask: torch.Tensor,
        reconstruction_indices: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        """Get reconstruction loss."""
        fill_in_the_blanks = self.fill_in_the_blanks_ratio > 0.0 and self.training

        if self.reconstruction_strategy == "dense" or reconstruction_indices is None:
            pass
        elif self.reconstruction_strategy.startswith("random_batch@"):
            reconstruction_indices = reconstruction_indices.to(x.device)
            x = x[:, reconstruction_indices]
            if fill_in_the_blanks:
                x_mask = x_mask[:, reconstruction_indices]
        else:
            raise NotImplementedError(
                f"Reconstruction strategy {self.reconstruction_strategy} not implemented."
            )

        if fill_in_the_blanks:
            reconst_loss = -(px.log_prob(x) * (1 - x_mask)).sum(dim=-1)
            mse = (
                torch.nn.functional.mse_loss(x * x_mask, px.mean * x_mask, reduction="none")
                .sum(dim=1)
                .mean(dim=0)
            )
        else:
            reconst_loss = -px.log_prob(x).sum(dim=-1)
            mse = torch.nn.functional.mse_loss(x, px.mean, reduction="none").sum(dim=1).mean(dim=0)

        return {
            MODULE_KEYS.RECONSTRUCTION_LOSS_KEY: reconst_loss,
            MODULE_KEYS.MSE_LOSS_KEY: mse,
        }

    def _get_kl_divergence_z(self, qz_m: torch.Tensor, qz_v: torch.Tensor) -> torch.Tensor:
        """Get KL divergence term for z."""
        qz = Normal(qz_m, torch.sqrt(qz_v))
        pz = Normal(torch.zeros_like(qz.mean), torch.ones_like(qz.mean))
        if self.prior == "normal":
            return kl_divergence(qz, pz).sum(dim=-1)
        else:
            raise NotImplementedError()

    def loss(
        self,
        tensors: TensorDict,
        inference_outputs: dict[str, Any],
        generative_outputs: dict[str, Any],
        kl_weight: float = 1.0,
    ) -> LossOutput:
        """Loss function.

        Parameters
        ----------
        tensors
            Dictionary containing tensor data.
        inference_outputs
            Outputs from the inference step.
        generative_outputs
            Outputs from the generative step.
        kl_weight
            Weight for KL divergence term.

        Returns
        -------
        LossOutput
            Loss output object containing various loss components.
        """
        x = tensors[REGISTRY_KEYS.X_KEY]
        x_mask = inference_outputs[MODULE_KEYS.X_MASK_KEY]
        qz_m = inference_outputs[MODULE_KEYS.QZM_KEY]
        qz_v = inference_outputs[MODULE_KEYS.QZV_KEY]
        px = generative_outputs[MODULE_KEYS.PX_KEY]
        reconstruction_indices = generative_outputs[MODULE_KEYS.RECONSTRUCTION_INDICES]

        kl_divergence_z = self._get_kl_divergence_z(qz_m, qz_v)
        reconst_losses = self._get_reconstruction_loss(px, x, x_mask, reconstruction_indices)
        reconst_loss = reconst_losses[MODULE_KEYS.RECONSTRUCTION_LOSS_KEY]
        assert kl_divergence_z.shape == reconst_loss.shape

        kl_local_for_warmup = kl_divergence_z
        kl_local_no_warmup = 0.0

        weighted_kl_local = kl_weight * kl_local_for_warmup + kl_local_no_warmup

        loss = torch.mean(reconst_loss + weighted_kl_local)

        kl_local = {MODULE_KEYS.KL_Z_KEY: kl_divergence_z}
        reconstruction_loss = {MODULE_KEYS.RECONSTRUCTION_LOSS_KEY: reconst_loss}
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            kl_local=kl_local,
            extra_metrics={
                MODULE_KEYS.MSE_LOSS_KEY: reconst_losses[MODULE_KEYS.MSE_LOSS_KEY],
            },
        )

    @torch.no_grad()
    def sample(
        self,
        tensors: TensorDict,
        n_samples: int = 1,
        library_size: int = 1,
        generative_kwargs: dict | None = None,
    ) -> torch.Tensor:
        # Note: Not tested
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Dictionary containing tensor data.
        n_samples
            Number of required samples for each cell.
        library_size
            Library size to scale samples to.
        generative_kwargs
            Keyword args for ``generative()`` in fwd pass

        Returns
        -------
        torch.Tensor
            Tensor with shape (n_cells, n_genes, n_samples).
        """
        inference_kwargs = dict(n_samples=n_samples)  # noqa: C408
        (
            _,
            generative_outputs,
        ) = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            generative_kwargs=generative_kwargs,
            compute_loss=False,
        )

        dist = generative_outputs[MODULE_KEYS.PX_KEY]

        if n_samples > 1:
            exprs = dist.sample().movedim(0, -1)
        else:
            exprs = dist.sample()

        return exprs.cpu()

    @torch.no_grad()
    @auto_move_data
    def marginal_ll(self, tensors: TensorDict, n_mc_samples: int) -> float:
        """Compute marginal log-likelihood.

        Parameters
        ----------
        tensors
            Dictionary containing tensor data.
        n_mc_samples
            Number of Monte Carlo samples for estimation.

        Returns
        -------
        float
            Marginal log-likelihood value.
        """
        sample_batch = tensors[REGISTRY_KEYS.X_KEY]

        to_sum = torch.zeros(sample_batch.size()[0], n_mc_samples)

        for i in range(n_mc_samples):
            # Distribution parameters and sampled variables
            inference_outputs, _, losses = self.forward(tensors)
            qz_m = inference_outputs[MODULE_KEYS.QZM_KEY]
            qz_v = inference_outputs[MODULE_KEYS.QZV_KEY]
            z = inference_outputs[MODULE_KEYS.Z_KEY]

            # Reconstruction Loss
            reconst_loss = losses.dict_sum(losses.reconstruction_loss)

            # Log-probabilities

            p_z = Normal(torch.zeros_like(qz_m), torch.ones_like(qz_v)).log_prob(z).sum(dim=-1)
            p_x_zl = -reconst_loss

            to_sum[:, i] = p_z + p_x_zl

        batch_log_lkl = torch.logsumexp(to_sum, dim=-1) - np.log(n_mc_samples)
        log_lkl = torch.sum(batch_log_lkl).item()
        return log_lkl
