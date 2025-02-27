from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import BaseModuleClass, EmbeddingModuleMixin, LossOutput, auto_move_data

from ._base_components import EncoderDecoder
from ._priors import StandardPrior, VampPrior

if TYPE_CHECKING:
    from typing import Literal

    from torch.distributions import Distribution

torch.backends.cudnn.benchmark = True


class SysVAE(BaseModuleClass, EmbeddingModuleMixin):
    """CVAE with optional VampPrior and latent cycle consistency loss.

    Described in
    `Hrovatin et al. (2023) <https://doi.org/10.1101/2023.11.03.565463>`_.

    Parameters
    ----------
    n_input
        Number of input features.
    n_batch
        Number of batches.
    n_continuous_cov
        Number of continuous covariates.
    n_cats_per_cov
        A list of integers containing the number of categories
        for each categorical covariate.
    embed_categorical_covariates
        If ``True`` embeds categorical covariates and batches
        into continuously-valued vectors instead of using one-hot encoding.
    prior
        Which prior distribution to use.
        * ``'standard_normal'``: Standard normal distribution.
        * ``'vamp'``: VampPrior.
    n_prior_components
        Number of prior components for VampPrior.
    trainable_priors
        Should prior components of VampPrior be trainable.
    pseudoinput_data
        Initialisation data for VampPrior.
        Should match input tensors structure.
    n_latent
        Numer of latent space dimensions.
    n_hidden
        Numer of hidden nodes per layer for encoder and decoder.
    n_layers
        Number of hidden layers for encoder and decoder.
    dropout_rate
        Dropout rate for encoder and decoder.
    out_var_mode
        How variance is predicted in decoder,
        see :class:`~scvi.external.sysvi.nn.VarEncoder`.
        One of the following:
        * ``'sample_feature'`` - learn variance per sample and feature.
        * ``'feature'`` - learn variance per feature, constant across samples.
    encoder_decoder_kwargs
        Additional kwargs passed to encoder and decoder.
        Both encoder and decoder use :class:`~scvi.external.sysvi.nn.EncoderDecoder`.
    embedding_kwargs
        Keyword arguments passed into :class:`~scvi.nn.Embedding`
        if ``embed_categorical_covariates`` is set to ``True``.
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int,
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        embed_categorical_covariates: bool = False,
        prior: Literal["standard_normal", "vamp"] = "vamp",
        n_prior_components: int = 5,
        trainable_priors: bool = True,
        pseudoinput_data: dict[str, torch.Tensor] | None = None,
        n_latent: int = 15,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        out_var_mode: Literal["sample_feature", "feature"] = "feature",
        encoder_decoder_kwargs: dict | None = None,
        embedding_kwargs: dict | None = None,
    ):
        super().__init__()

        self.embed_categorical_covariates = embed_categorical_covariates

        encoder_decoder_kwargs = encoder_decoder_kwargs or {}
        embedding_kwargs = embedding_kwargs or {}

        self.n_batch = n_batch
        n_cat_list = [n_batch]
        n_cont = n_continuous_cov
        if n_cats_per_cov is not None:
            if self.embed_categorical_covariates:
                for idx, n in enumerate(n_cats_per_cov):
                    covariate_name = f"cov{idx}"
                    self.init_embedding(covariate_name, n, **embedding_kwargs)
                    n_cont += self.get_embedding(covariate_name).embedding_dim
            else:
                n_cat_list.extend(n_cats_per_cov)

        self.encoder = EncoderDecoder(
            n_input=n_input,
            n_output=n_latent,
            n_cat_list=n_cat_list,
            n_cont=n_cont,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            sample=True,
            var_mode="sample_feature",
            **encoder_decoder_kwargs,
        )

        self.decoder = EncoderDecoder(
            n_input=n_latent,
            n_output=n_input,
            n_cat_list=n_cat_list,
            n_cont=n_cont,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            sample=True,
            var_mode=out_var_mode,
            **encoder_decoder_kwargs,
        )

        if prior == "standard_normal":
            self.prior = StandardPrior()
        elif prior == "vamp":
            assert pseudoinput_data is not None, (
                "Pseudoinput data must be specified if using VampPrior"
            )
            pseudoinput_data = self._get_inference_input(pseudoinput_data)
            self.prior = VampPrior(
                n_components=n_prior_components,
                encoder=self.encoder,
                x=pseudoinput_data["x"],
                n_cat_list=n_cat_list,
                batch_index=pseudoinput_data["batch_index"],
                cat_list=pseudoinput_data["cat_covs"],
                cont=pseudoinput_data["cont_covs"],
                trainable_priors=trainable_priors,
            )
        else:
            raise ValueError("Prior not recognised")

    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor], **kwargs
    ) -> dict[str, torch.Tensor | list[torch.Tensor] | None]:
        """Parse the input tensors to get inference inputs.

        Parameters
        ----------
        tensors
            Input tensors.
        kwargs
            Not used. Added for inheritance compatibility.
        """
        cov = self._get_cov(tensors=tensors)
        return {
            MODULE_KEYS.X_KEY: tensors[REGISTRY_KEYS.X_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            MODULE_KEYS.CONT_COVS_KEY: cov["continuous"],
            MODULE_KEYS.CAT_COVS_KEY: cov["categorical"],
        }

    def _get_inference_cycle_input(
        self,
        tensors: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        cycle_batch: torch.Tensor,
        **kwargs,
    ) -> dict[str, torch.Tensor | list[torch.Tensor] | None]:
        """Input tensors, generative outputs, and cycle batch -> cycle inference inputs.

        Parameters
        ----------
        tensors
            Input tensors.
        generative_outputs
            Outputs of the generative pass.
        cycle_batch
            Batch covariate to be used for the cycle inference.
            dim = n_samples * 1
        kwargs
            Not used. Added for inheritance compatibility.

        Note: cycle covariates differ from the original publication.
        Instead of mock covariates the real input covariates are used in cycle.

        """
        cov = self._get_cov(tensors=tensors)
        return {
            MODULE_KEYS.X_KEY: generative_outputs["px"].loc,
            MODULE_KEYS.BATCH_INDEX_KEY: cycle_batch,
            MODULE_KEYS.CONT_COVS_KEY: cov["continuous"],
            MODULE_KEYS.CAT_COVS_KEY: cov["categorical"],
        }

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        cycle_batch: torch.Tensor,
        **kwargs,
    ) -> dict[str, torch.Tensor | dict[str, torch.Tensor | list[torch.Tensor] | None]]:
        """Input tensors, inference outputs, and cycle batch info -> generative inputs.

        Parameters
        ----------
        tensors
            Input tensors.
        inference_outputs
            Outputs of the inference pass.
        cycle_batch
            Batch covariate to be used for the cycle expression generation.
            dim = n_samples * 1
        kwargs
            Not used. Added for inheritance compatibility.

        Note: cycle covariates differ from the original publication.
        Instead of mock covariates the real input covariates are used in cycle.
        """
        cov = self._get_cov(tensors=tensors)
        return {
            MODULE_KEYS.Z_KEY: inference_outputs[MODULE_KEYS.Z_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            "cycle_batch": cycle_batch,
            MODULE_KEYS.CONT_COVS_KEY: cov["continuous"],
            MODULE_KEYS.CAT_COVS_KEY: cov["categorical"],
        }

    @auto_move_data
    def _get_cov(
        self,
        tensors: dict[str, torch.Tensor],
    ) -> dict[str, torch.Tensor | list[torch.Tensor] | None]:
        """Process all covs into continuous and categorical components for cVAE.

        Parameters
        ----------
        tensors
            Input tensors.

        Returns
        -------
        Covariates that can be used for decoder and encoder.
        Keys:
        * ``'cat'``: All covariates that require one-hot encoding.
                     List of tensors with dim = n_samples * 1.
                     If absent returns empty list.
        * ``'cont'``: All covariates that are already continous.
                      Includes continous and embedded categorical covariates.
                      Single tensor of dim = n_samples * n_concatenated_cov_features.
                      If absent returns None.
        """
        cat_parts = []
        cont_parts = []
        if REGISTRY_KEYS.CONT_COVS_KEY in tensors:
            cont_parts.append(tensors[REGISTRY_KEYS.CONT_COVS_KEY])
        if REGISTRY_KEYS.CAT_COVS_KEY in tensors:
            cat = torch.split(tensors[REGISTRY_KEYS.CAT_COVS_KEY], 1, dim=1)
            if self.embed_categorical_covariates:
                for idx, tensor in enumerate(cat):
                    cont_parts.append(self.compute_embedding(f"cov{idx}", tensor))
            else:
                cat_parts.extend(cat)
        cov = {
            "categorical": cat_parts,
            "continuous": torch.concat(cont_parts, dim=1) if len(cont_parts) > 0 else None,
        }
        return cov

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        batch_index: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, torch.Tensor | Distribution | None]:
        """Inference: expression & covariates -> latent representation."""
        result = self.encoder(x, batch_index=batch_index, cat_list=cat_covs, cont=cont_covs)
        z, qz = result["q"], result["q_dist"]
        return {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.QZ_KEY: qz,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        batch_index: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: list[torch.Tensor] | None = None,
        cycle_batch: torch.Tensor | None = None,
        compute_original: bool = True,
        compute_cycle: bool = True,
        y: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        """Generative: latent representation & covariates -> expression."""
        res = {}
        if compute_original:
            res["px"] = self.decoder(
                x=z, batch_index=batch_index, cont=cont_covs, cat_list=cat_covs
            )["q_dist"]
        if compute_cycle:
            res["px_cycle"] = self.decoder(
                x=z, batch_index=cycle_batch, cont=cont_covs, cat_list=cat_covs
            )["q_dist"]
        return res

    @auto_move_data
    def forward(
        self,
        tensors: dict[str, torch.Tensor],
        get_inference_input_kwargs: dict | None = None,
        get_generative_input_kwargs: dict | None = None,
        inference_kwargs: dict | None = None,
        generative_kwargs: dict | None = None,
        loss_kwargs: dict | None = None,
        compute_loss: bool = True,
    ) -> (
        tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]
        | tuple[dict[str, torch.Tensor], dict[str, torch.Tensor], LossOutput]
    ):
        # Parse kwargs
        inference_kwargs = inference_kwargs or {}
        generative_kwargs = generative_kwargs or {}
        loss_kwargs = loss_kwargs or {}
        get_inference_input_kwargs = get_inference_input_kwargs or {}
        get_generative_input_kwargs = get_generative_input_kwargs or {}

        # Inference
        inference_inputs = self._get_inference_input(tensors, **get_inference_input_kwargs)
        inference_outputs = self.inference(**inference_inputs, **inference_kwargs)
        # Generative
        cycle_batch = self.random_select_batch(tensors[REGISTRY_KEYS.BATCH_KEY])
        generative_inputs = self._get_generative_input(
            tensors,
            inference_outputs,
            cycle_batch=cycle_batch,
            **get_generative_input_kwargs,
        )
        if loss_kwargs.get("z_distance_cycle_weight", 2) == 0:
            compute_cycle = False
        else:
            compute_cycle = True
        generative_outputs = self.generative(
            **generative_inputs,
            compute_original=True,
            compute_cycle=compute_cycle,
            **generative_kwargs,
        )
        if compute_cycle:
            # Inference cycle
            inference_cycle_inputs = self._get_inference_cycle_input(
                tensors=tensors,
                generative_outputs=generative_outputs,
                cycle_batch=cycle_batch,
                **get_inference_input_kwargs,
            )
            inference_cycle_outputs = self.inference(**inference_cycle_inputs, **inference_kwargs)
            inference_outputs["qz_cycle"] = inference_cycle_outputs["qz"]

        if compute_loss:
            losses = self.loss(
                tensors=tensors,
                inference_outputs=inference_outputs,
                generative_outputs=generative_outputs,
                compute_cycle=compute_cycle,
                **loss_kwargs,
            )
            return inference_outputs, generative_outputs, losses
        else:
            return inference_outputs, generative_outputs

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        kl_weight: float = 1.0,
        reconstruction_weight: float = 1.0,
        z_distance_cycle_weight: float = 2.0,
        compute_cycle: bool = True,
    ) -> LossOutput:
        """Compute loss of forward pass.

        Parameters
        ----------
        tensors
            Input tensors.
        inference_outputs
            Outputs of normal and cycle inference pass.
        generative_outputs
            Outputs of the normal generative pass.
        kl_weight
            Weight for KL loss.
        reconstruction_weight
            Weight for reconstruction loss.
        z_distance_cycle_weight
            Weight for cycle loss.

        Returns
        -------
        Loss components:
        Cycle loss is added to extra metrics as ``'cycle_loss'``.
        """
        # Reconstruction loss
        x_true = tensors[REGISTRY_KEYS.X_KEY]
        reconst_loss_x = torch.nn.GaussianNLLLoss(reduction="none")(
            generative_outputs["px"].loc, x_true, generative_outputs["px"].scale
        ).sum(dim=1)
        reconst_loss = reconst_loss_x

        # Kl divergence on latent space
        kl_divergence_z = self.prior.kl(
            qz=inference_outputs["qz"],
            z=inference_outputs["z"],
        )

        if z_distance_cycle_weight > 0:
            z_distance_cyc = self.latent_cycle_consistency(
                qz=inference_outputs["qz"], qz_cycle=inference_outputs["qz_cycle"]
            )
            if "batch_weights" in tensors.keys():
                z_distance_cyc *= tensors["batch_weights"].flatten()
        else:
            z_distance_cyc = torch.zeros_like(reconst_loss)

        loss = (
            reconst_loss * reconstruction_weight
            + kl_divergence_z * kl_weight
            + z_distance_cyc * z_distance_cycle_weight
        )

        return LossOutput(
            loss=loss.mean(),
            reconstruction_loss=reconst_loss,
            kl_local=kl_divergence_z,
            extra_metrics={"cycle_loss": z_distance_cyc.mean()},
        )

    def random_select_batch(self, batch: torch.Tensor) -> torch.Tensor:
        """Randomly selects a new batch different from the real one for each cell.

        Parameters
        ----------
        batch : torch.Tensor
            Tensor containing the real batch index for each cell.

        Returns
        -------
        torch.Tensor
            Tensor with newly assigned batch indices for each cell.
        """
        batch_one_hot = torch.nn.functional.one_hot(batch.squeeze(-1), self.n_batch)
        available_batches = torch.nonzero(1 - batch_one_hot, as_tuple=True)[1]
        available_batches = available_batches.view(batch.shape[0], -1)

        # Select a new batch index for each cell from the available batches
        new_batch_indices = available_batches[
            torch.arange(batch.shape[0], device=batch.device),
            torch.randint(0, available_batches.shape[1], (batch.shape[0],), device=batch.device),
        ]

        return new_batch_indices.unsqueeze(-1)

    @staticmethod
    def latent_cycle_consistency(
        qz: torch.Tensor,
        qz_cycle: torch.Tensor,
    ) -> torch.Tensor:
        """MSE loss between standardised inputs.

        MSE loss should be computed on standardized latent representations
        as else model can learn to cheat the MSE loss
        by setting the latent representations to smaller numbers.
        Standardizer is fitted on concatenation of both inputs.

        Parameters
        ----------
        qz
            Posterior distribution from the inference pass.
        qz_cycle
            Posterior distribution from the cycle inference pass.
        """
        # Standardise data (jointly both z-s) before MSE calculation
        z = torch.concat([qz.loc, qz_cycle.loc])
        means = z.mean(dim=0, keepdim=True)
        stds = z.std(dim=0, keepdim=True)

        def standardize(x: torch.Tensor) -> torch.Tensor:
            return (x - means) / stds

        return torch.nn.MSELoss(reduction="none")(
            standardize(qz.loc), standardize(qz_cycle.loc)
        ).sum(dim=1)

    @torch.inference_mode()
    def sample(self, *args, **kwargs):
        """Generate expression samples from posterior generative distribution.

        Not implemented as the use of decoded expression
        is not recommended for SysVI.

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError("The use of decoded expression is not recommended for SysVI.")
