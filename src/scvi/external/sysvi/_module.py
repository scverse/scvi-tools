from __future__ import annotations

import torch
from typing_extensions import Literal

from scvi import REGISTRY_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

from ._base_components import Embedding, EncoderDecoder
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
    n_system
        Number of systems.
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
        n_cov_const: int,
        cov_embed_sizes: list,
        n_system: int,
        cov_embed_dims: int = 10,
        prior: Literal["standard_normal", "vamp"] = "vamp",
        n_prior_components: int = 5,
        trainable_priors: bool = True,
        pseudoinput_data: dict[str, torch.Tensor] | None = None,
        n_latent: int = 15,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        out_var_mode: str = "feature",
        **enc_dec_kwargs,
    ):
        super().__init__()

        self.embed_cov = len(cov_embed_sizes) > 0  # Will any covs be embedded
        self.n_cov_const = n_cov_const  # Dimension of covariates that are not embedded
        n_cov = (
            n_cov_const + len(cov_embed_sizes) * cov_embed_dims
        )  # Total size of covs (embedded & not embedded)
        n_cov_encoder = n_cov + n_system  # N covariates passed to Module (cov & system)

        if self.embed_cov:
            self.cov_embeddings = torch.nn.ModuleList(
                [Embedding(size=size, cov_embed_dims=cov_embed_dims) for size in cov_embed_sizes]
            )

        self.encoder = EncoderDecoder(
            n_input=n_input,
            n_output=n_latent,
            n_cov=n_cov_encoder,
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
            n_cov=n_cov_encoder,
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
            if pseudoinput_data is not None:
                pseudoinput_data = self._get_inference_input(pseudoinput_data)
            self.prior = VampPrior(
                n_components=n_prior_components,
                n_input=n_input,
                n_cov=n_cov_encoder,
                encoder=self.encoder,
                data=(
                    pseudoinput_data["expr"],
                    self._merge_cov(
                        cov=pseudoinput_data["cov"], system=pseudoinput_data["system"]
                    ),
                ),
                trainable_priors=trainable_priors,
            )
        else:
            raise ValueError("Prior not recognised")

    def _get_inference_input(self, tensors, **kwargs) -> dict[str, torch.Tensor]:
        """Parse the input tensors to get inference inputs"""
        expr = tensors[REGISTRY_KEYS.X_KEY]
        cov = self._get_cov(tensors=tensors)
        system = tensors["system"]
        input_dict = {"expr": expr, "cov": cov, "system": system}
        return input_dict

    def _get_inference_cycle_input(
        self, tensors, generative_outputs, selected_system: torch.Tensor, **kwargs
    ) -> dict[str, torch.Tensor]:
        """Parse the input tensors and cycle system info to get cycle inference inputs"""
        expr = generative_outputs["y_m"]
        cov = self._mock_cov(self._get_cov(tensors=tensors))
        system = selected_system
        input_dict = {"expr": expr, "cov": cov, "system": system}
        return input_dict

    def _get_generative_input(
        self, tensors, inference_outputs, selected_system: torch.Tensor, **kwargs
    ) -> dict[str, torch.Tensor]:
        """Parse the input tensors, inference inputs, and cycle system to get generative inputs"""
        z = inference_outputs["z"]

        cov = self._get_cov(tensors=tensors)
        cov = {"x": cov, "y": self._mock_cov(cov)}

        system = {"x": tensors["system"], "y": selected_system}

        input_dict = {"z": z, "cov": cov, "system": system}
        return input_dict

    @auto_move_data
    def _get_cov(self, tensors: dict[str, torch.Tensor]) -> torch.Tensor | None:
        """Merge all covariates into single tensor, including embedding of covariates"""
        cov = []
        if self.n_cov_const > 0:
            cov.append(tensors["covariates"])
        if self.embed_cov:
            cov.extend(
                [
                    embedding(tensors["covariates_embed"][:, idx].int())
                    for idx, embedding in enumerate(self.cov_embeddings)
                ]
            )
        cov = torch.concat(cov, dim=1) if len(cov) > 0 else None
        return cov

    @staticmethod
    def _merge_cov(cov: torch.Tensor | None, system: torch.Tensor) -> torch.Tensor:
        """Merge full covariate data and system data to get cov for model input"""
        return torch.cat([cov, system], dim=1) if cov is not None else system

    @staticmethod
    def _mock_cov(cov: torch.Tensor | None) -> torch.Tensor | None:
        """Make mock (all 0) covariates for cycle"""
        return torch.zeros_like(cov) if cov is not None else None

    @auto_move_data
    def inference(self, expr, cov, system) -> dict:
        """
        expression & cov -> latent representation

        Parameters
        ----------
        expr
            Expression data
        cov
            Full covariate data (categorical, categorical embedded, and continuous
        system
            System representation

        Returns
        -------
        Posterior parameters and sample
        """
        z = self.encoder(x=expr, cov=self._merge_cov(cov=cov, system=system))
        return {"z": z["y"], "z_m": z["y_m"], "z_v": z["y_v"]}

    @auto_move_data
    def generative(self, z, cov, system, x_x: bool = True, x_y: bool = True) -> dict:
        """
        latent representation & cov -> expression

        Parameters
        ----------
        z
            Latent embedding
        cov
            Full covariate data (categorical, categorical embedded, and continuous
        system
            System representation
        x_x
            Decode to original system
        x_y
            Decode to replacement system

        Returns
        -------
        Decoded distribution parameters and sample
        """

        def outputs(compute, name, res, x, cov, system):
            if compute:
                res_sub = self.decoder(x=x, cov=self._merge_cov(cov=cov, system=system))
                res[name] = res_sub["y"]
                res[name + "_m"] = res_sub["y_m"]
                res[name + "_v"] = res_sub["y_v"]

        res = {}
        outputs(compute=x_x, name="x", res=res, x=z, cov=cov["x"], system=system["x"])
        outputs(compute=x_y, name="y", res=res, x=z, cov=cov["y"], system=system["y"])
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
        inference_kwargs = _get_dict_if_none(inference_kwargs)
        generative_kwargs = _get_dict_if_none(generative_kwargs)
        loss_kwargs = _get_dict_if_none(loss_kwargs)
        get_inference_input_kwargs = _get_dict_if_none(get_inference_input_kwargs)
        get_generative_input_kwargs = _get_dict_if_none(get_generative_input_kwargs)

        # Inference
        inference_inputs = self._get_inference_input(tensors, **get_inference_input_kwargs)
        inference_outputs = self.inference(**inference_inputs, **inference_kwargs)
        # Generative
        selected_system = self.random_select_systems(tensors["system"])
        generative_inputs = self._get_generative_input(
            tensors,
            inference_outputs,
            selected_system=selected_system,
            **get_generative_input_kwargs,
        )
        generative_outputs = self.generative(
            **generative_inputs, x_x=True, x_y=True, **generative_kwargs
        )
        # Inference cycle
        inference_cycle_inputs = self._get_inference_cycle_input(
            tensors=tensors,
            generative_outputs=generative_outputs,
            selected_system=selected_system,
            **get_inference_input_kwargs,
        )
        inference_cycle_outputs = self.inference(**inference_cycle_inputs, **inference_kwargs)

        # Combine outputs of all forward pass components
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
    ):
        x_true = tensors[REGISTRY_KEYS.X_KEY]

        # Reconstruction loss
        reconst_loss_x = torch.nn.GaussianNLLLoss(reduction="none")(
            generative_outputs["x_m"], x_true, generative_outputs["x_v"]).sum(dim=1)

        reconst_loss = reconst_loss_x

        # Kl divergence on latent space
        kl_divergence_z = self.prior.kl(
            m_q=inference_outputs["z_m"],
            v_q=inference_outputs["z_v"],
            z=inference_outputs["z"],
        )

        def z_dist(z_x_m: torch.Tensor, z_y_m: torch.Tensor):
            """MSE loss between standardised inputs with standardizer fitted on concatenation of both inputs"""
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

        loss = (
            reconst_loss * reconstruction_weight
            + kl_divergence_z * kl_weight
            + z_distance_cyc * z_distance_cycle_weight
        )

        return LossOutput(
            n_obs_minibatch=loss.shape[0],
            loss=loss.mean(),
            extra_metrics={
                "reconstruction_loss": reconst_loss.mean(),
                "kl_local": kl_divergence_z.mean(),
                "z_distance_cycle": z_distance_cyc.mean(),
            },
        )

    @staticmethod
    def random_select_systems(system: torch.Tensor) -> torch.Tensor:
        """For every cell randomly selects a new system that is different from the original system

        Parameters
        ----------
        system
            One hot encoded system information for each cell

        Returns
        -------
        One hot encoding of newly selected system for each cell

        """
        # get available systems - those that are zero will become nonzero and vice versa
        available_systems = 1 - system
        # Get nonzero indices for each cell
        row_indices, col_indices = torch.nonzero(available_systems, as_tuple=True)
        col_pairs = col_indices.view(-1, system.shape[1] - 1)
        # Select system for every cell from available systems
        randomly_selected_indices = col_pairs.gather(
            1,
            torch.randint(
                0,
                system.shape[1] - 1,
                size=(col_pairs.size(0), 1),
                device=col_pairs.device,
                dtype=col_pairs.dtype,
            ),
        )
        new_tensor = torch.zeros_like(available_systems)
        # generate system covariate tensor
        new_tensor.scatter_(1, randomly_selected_indices, 1)

        return new_tensor

    def sample(self, *args, **kwargs):
        raise NotImplementedError("")


def _get_dict_if_none(param: dict | None) -> dict:
    """If not a dict return empty dict"""
    param = {} if not isinstance(param, dict) else param
    return param
