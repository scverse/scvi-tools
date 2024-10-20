from __future__ import annotations

from typing import Literal

import torch

from scvi import REGISTRY_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

from ._base_components import EncoderDecoder
from ._priors import StandardPrior, VampPrior

torch.backends.cudnn.benchmark = True


class SysVAE(BaseModuleClass):
    """CVAE with optional VampPrior and latent cycle consistency loss.

    Parameters
    ----------
    n_input
        Number of input features.
        Passed directly from Model.
    n_cov_const
        Dimensionality of covariate data that will not be further embedded.
        Passed directly from Model.
    cov_embed_sizes
        Number of categories per every cov to be embedded, e.g. [cov1_n_categ, cov2_n_categ, ...].
        Passed directly from Model.
    n_batch
        Number of batches.
        Passed directly from Model.
    cov_embed_dims
        Dimension for covariate embedding.
    prior
        Which prior distribution to use.
        Passed directly from Model.
    n_prior_components
        If VampPrior - how many prior components to use.
        Passed directly from Model.
    trainable_priors
        If VampPrior - should prior components be trainable.
    pseudoinput_data
        Initialisation data for VampPrior. Should match input tensors structure.
        Passed directly from Model.
    n_latent
        Numer of latent space dimensions.
    n_hidden
        Number of nodes in hidden layers.
    n_layers
        Number of hidden layers.
    dropout_rate
        Dropout rate.
    out_var_mode
        See :class:`~scvi.external.sysvi.nn.VarEncoder`
    enc_dec_kwargs
        Additional kwargs passed to encoder and decoder.
    """

    # TODO could disable computation of cycle if predefined that cycle wil not be used

    def __init__(
        self,
        n_input: int,
        n_batch: int,
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        embed_cat: bool = False,
        prior: Literal["standard_normal", "vamp"] = "vamp",
        n_prior_components: int = 5,
        trainable_priors: bool = True,
        pseudoinput_data: dict[str, torch.Tensor] | None = None,
        n_latent: int = 15,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        out_var_mode: Literal["sample_feature", "feature"] = "feature",
        enc_dec_kwargs: dict | None = None,
        embedding_kwargs: dict | None = None,
    ):
        super().__init__()

        self.embed_cat = embed_cat

        enc_dec_kwargs = enc_dec_kwargs or {}
        embedding_kwargs = embedding_kwargs or {}

        self.n_batch = n_batch
        n_cat_list = [n_batch]
        n_cont = n_continuous_cov
        if n_cats_per_cov is not None:
            if self.embed_cat:
                for cov, n in enumerate(n_cats_per_cov):
                    cov = self._cov_idx_name(cov=cov)
                    self.init_embedding(cov, n, **embedding_kwargs)
                    n_cont += self.get_embedding(cov).embedding_dim
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
            **enc_dec_kwargs,
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
            **enc_dec_kwargs,
        )

        if prior == "standard_normal":
            self.prior = StandardPrior()
        elif prior == "vamp":
            assert (
                pseudoinput_data is not None
            ), "Pseudoinput data must be specified if using VampPrior"
            pseudoinput_data = self._get_inference_input(pseudoinput_data)
            self.prior = VampPrior(
                n_components=n_prior_components,
                encoder=self.encoder,
                data_x=pseudoinput_data["expr"],
                n_cat_list=n_cat_list,
                data_cat=self._merge_batch_cov(
                    cat=pseudoinput_data["cat"], batch=pseudoinput_data["batch"]
                ),
                data_cont=pseudoinput_data["cont"],
                trainable_priors=trainable_priors,
            )
        else:
            raise ValueError("Prior not recognised")

    @staticmethod
    def _cov_idx_name(cov: int | float):
        return "cov" + str(cov)

    def _get_inference_input(self, tensors, **kwargs) -> dict[str, torch.Tensor]:
        """Parse the input tensors to get inference inputs"""
        cov = self._get_cov(tensors=tensors)
        input_dict = {
            "expr": tensors[REGISTRY_KEYS.X_KEY],
            "batch": tensors[REGISTRY_KEYS.BATCH_KEY],
            "cat": cov["cat"],
            "cont": cov["cont"],
        }
        return input_dict

    def _get_inference_cycle_input(
        self, tensors, generative_outputs, selected_batch: torch.Tensor, **kwargs
    ) -> dict[str, torch.Tensor]:
        """Parse the input tensors and cycle batch info to get cycle inference inputs."""
        cov = self._mock_cov(self._get_cov(tensors=tensors))
        input_dict = {
            "expr": generative_outputs["y_m"],
            "batch": selected_batch,
            "cat": cov["cat"],
            "cont": cov["cont"],
        }
        return input_dict

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        selected_batch: torch.Tensor,
        **kwargs,
    ) -> dict[str, torch.Tensor | dict[str, torch.Tensor | list[torch.Tensor] | None]]:
        """Parse the input tensors, inference inputs, and cycle batch to get generative inputs"""
        z = inference_outputs["z"]

        cov = self._get_cov(tensors=tensors)
        cov_mock = self._mock_cov(cov)
        cat = {"x": cov["cat"], "y": cov_mock["cat"]}
        cont = {"x": cov["cont"], "y": cov_mock["cont"]}

        batch = {"x": tensors["batch"], "y": selected_batch}

        input_dict = {"z": z, "batch": batch, "cat": cat, "cont": cont}
        return input_dict

    @auto_move_data  # TODO remove?
    def _get_cov(
        self, tensors: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor | list[torch.Tensor] | None]:
        """Process all covariates into continuous and categorical components for cVAE"""
        cat_parts = []
        cont_parts = []
        if REGISTRY_KEYS.CONT_COVS_KEY in tensors:
            cont_parts.append(tensors[REGISTRY_KEYS.CONT_COVS_KEY])
        if REGISTRY_KEYS.CAT_COVS_KEY in tensors:
            cat = torch.split(tensors[REGISTRY_KEYS.CAT_COVS_KEY], 1, dim=1)
            if self.embed_cat:
                for idx, tensor in enumerate(cat):
                    cont_parts.append(self.compute_embedding(tensor, self._cov_idx_name(idx)))
            else:
                cat_parts.extend(cat)
        cov = {
            "cat": cat_parts,
            "cont": torch.concat(cont_parts, dim=1) if len(cont_parts) > 0 else None,
        }
        return cov

    @staticmethod
    def _merge_batch_cov(cat: list[torch.Tensor], batch: torch.Tensor) -> list[torch.Tensor]:
        return [batch] + cat

    @staticmethod
    def _mock_cov(cov: dict[str, list[torch.Tensor], torch.Tensor, None]) -> torch.Tensor | None:
        """Make mock (all 0) covariates for cycle"""
        mock = {
            "cat": [torch.zeros_like(cat) for cat in cov["cat"]],
            "cont": torch.zeros_like(cov["cont"]) if cov["cont"] is not None else None,
        }
        return mock

    @auto_move_data
    def inference(
        self,
        expr: torch.Tensor,
        batch: torch.Tensor,
        cat: list[torch.Tensor],
        cont: torch.Tensor | None,
    ) -> dict:
        """
        expression & cov -> latent representation

        Parameters
        ----------
        expr
            Expression data
        cov
            Full covariate data (categorical, categorical embedded, and continuous
        batch
            System representation

        Returns
        -------
        Posterior parameters and sample
        """
        z = self.encoder(x=expr, cat_list=self._merge_batch_cov(cat=cat, batch=batch), cont=cont)
        return {"z": z["y"], "z_m": z["y_m"], "z_v": z["y_v"]}

    @auto_move_data
    def generative(self, z, batch, cat, cont, x_x: bool = True, x_y: bool = True) -> dict:
        """
        latent representation & cov -> expression

        Parameters
        ----------
        z
            Latent embedding
        cov
            Full covariate data (categorical, categorical embedded, and continuous
        batch
            System representation
        x_x
            Decode to original batch
        x_y
            Decode to replacement batch

        Returns
        -------
        Decoded distribution parameters and sample
        """

        def outputs(
            compute: bool,
            name: str,
            res: dict,
            x: torch.Tensor,
            batch: torch.Tensor,
            cat: list[torch.Tensor],
            cont: torch.Tensor | None,
        ):
            if compute:
                res_sub = self.decoder(
                    x=x, cat_list=self._merge_batch_cov(cat=cat, batch=batch), cont=cont
                )
                res[name] = res_sub["y"]
                res[name + "_m"] = res_sub["y_m"]
                res[name + "_v"] = res_sub["y_v"]

        res = {}
        outputs(
            compute=x_x, name="x", res=res, x=z, batch=batch["x"], cat=cat["x"], cont=cont["x"]
        )
        outputs(
            compute=x_y, name="y", res=res, x=z, batch=batch["y"], cat=cat["y"], cont=cont["y"]
        )
        return res

    @auto_move_data
    def forward(
        self,
        tensors,
        get_inference_input_kwargs: dict | None = None,
        get_generative_input_kwargs: dict | None = None,
        inference_kwargs: dict | None = None,
        generative_kwargs: dict | None = None,
        loss_kwargs: dict | None = None,
        compute_loss=True,
    ) -> (
        tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]
        | tuple[dict[str, torch.Tensor], dict[str, torch.Tensor], LossOutput]
    ):
        """
        Forward pass through the network.

        Parameters
        ----------
        tensors
            tensors to pass through
        get_inference_input_kwargs
            Keyword args for ``_get_inference_input()``
        get_generative_input_kwargs
            Keyword args for ``_get_generative_input()``
        inference_kwargs
            Keyword args for ``inference()``
        generative_kwargs
            Keyword args for ``generative()``
        loss_kwargs
            Keyword args for ``loss()``
        compute_loss
            Whether to compute loss on forward pass. This adds
            another return value.
        """
        """Core of the forward call shared by PyTorch- and Jax-based modules."""

        # TODO currently some forward paths are computed despite potentially having loss weight=0 -
        #  don't compute if not needed

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
        selected_batch = self.random_select_batch(tensors[REGISTRY_KEYS.BATCH_KEY])
        generative_inputs = self._get_generative_input(
            tensors,
            inference_outputs,
            selected_batch=selected_batch,
            **get_generative_input_kwargs,
        )
        generative_outputs = self.generative(
            **generative_inputs, x_x=True, x_y=True, **generative_kwargs
        )
        # Inference cycle
        inference_cycle_inputs = self._get_inference_cycle_input(
            tensors=tensors,
            generative_outputs=generative_outputs,
            selected_batch=selected_batch,
            **get_inference_input_kwargs,
        )
        inference_cycle_outputs = self.inference(**inference_cycle_inputs, **inference_kwargs)

        # Combine outputs of all forward pass components - first and cycle pass
        inference_outputs_merged = dict(**inference_outputs)
        inference_outputs_merged.update(
            **{k.replace("z", "z_cyc"): v for k, v in inference_cycle_outputs.items()}
        )
        generative_outputs_merged = dict(**generative_outputs)

        if compute_loss:
            losses = self.loss(
                tensors=tensors,
                inference_outputs=inference_outputs_merged,
                generative_outputs=generative_outputs_merged,
                **loss_kwargs,
            )
            return inference_outputs_merged, generative_outputs_merged, losses
        else:
            return inference_outputs_merged, generative_outputs_merged

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        reconstruction_weight: float = 1.0,
        z_distance_cycle_weight: float = 2.0,
    ) -> LossOutput:
        # Reconstruction loss
        x_true = tensors[REGISTRY_KEYS.X_KEY]
        reconst_loss_x = torch.nn.GaussianNLLLoss(reduction="none")(
            generative_outputs["x_m"], x_true, generative_outputs["x_v"]
        ).sum(dim=1)
        reconst_loss = reconst_loss_x

        # Kl divergence on latent space
        kl_divergence_z = self.prior.kl(
            m_q=inference_outputs["z_m"],
            v_q=inference_outputs["z_v"],
            z=inference_outputs["z"],
        )

        def z_dist(z_x_m: torch.Tensor, z_y_m: torch.Tensor):
            """MSE loss between standardised inputs with standardizer fitted on concatenation of both inputs

            MSE loss should be computed on standardized latent values as else model can learn to cheat the MSE
            loss by putting latent parameters to even smaller numbers.
            """
            # Standardise data (jointly both z-s) before MSE calculation
            z = torch.concat([z_x_m, z_y_m])
            means = z.mean(dim=0, keepdim=True)
            stds = z.std(dim=0, keepdim=True)

            def standardize(x):
                return (x - means) / stds

            return torch.nn.MSELoss(reduction="none")(standardize(z_x_m), standardize(z_y_m)).sum(
                dim=1
            )

        z_distance_cyc = z_dist(z_x_m=inference_outputs["z_m"], z_y_m=inference_outputs["z_cyc_m"])
        if "batch_weights" in tensors.keys():
            z_distance_cyc *= tensors["batch_weights"].flatten()

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
        """For every cell randomly selects a new batch that is different from the original batch

        Parameters
        ----------
        batch
            One hot encoded batch information for each cell

        Returns
        -------
        One hot encoding of newly selected batch for each cell

        """
        # Get available batches -
        # those that are zero will become nonzero and vice versa
        batch = torch.nn.functional.one_hot(batch.squeeze(-1), self.n_batch)
        available_batches = 1 - batch
        # Get nonzero indices for each cell
        row_indices, col_indices = torch.nonzero(available_batches, as_tuple=True)
        col_pairs = col_indices.view(-1, batch.shape[1] - 1)
        # Select batch for every cell from available batches
        randomly_selected_indices = col_pairs.gather(
            1,
            torch.randint(
                0,
                batch.shape[1] - 1,
                size=(col_pairs.size(0), 1),
                device=col_pairs.device,
                dtype=col_pairs.dtype,
            ),
        )
        new_tensor = torch.zeros_like(available_batches)
        # generate batch covariate tensor
        new_tensor.scatter_(1, randomly_selected_indices, 1)

        return new_tensor

    @torch.inference_mode()
    def sample(self, *args, **kwargs):
        raise NotImplementedError("The use of decoded expression is not recommended for SysVI.")
