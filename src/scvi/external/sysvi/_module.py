from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

import torch

from scvi import REGISTRY_KEYS
from scvi.module.base import BaseModuleClass, EmbeddingModuleMixin, LossOutput, auto_move_data

from ._base_components import EncoderDecoder
from ._priors import StandardPrior, VampPrior

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
    embed_cat
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
    enc_dec_kwargs
        Additional kwargs passed to encoder and decoder.
        For both encoder and decoder
        :class:`~scvi.external.sysvi.nn.EncoderDecoder` is used.
    embedding_kwargs
        Keyword arguments passed into :class:`~scvi.nn.Embedding`
        if ``embed_cat`` is set to ``True``.
    """

    # TODO could disable computation of cycle if predefined
    #  that cycle loss will not be used.
    #  Cycle loss is not expected to be disabled in practice
    #  for typical use cases.
    #  As the use of cycle is currently only based on loss kwargs,
    #  which are specified only later, it can not be inferred here.

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
    def _cov_idx_name(cov: int) -> str:
        """Convert covariate index into a name used for embedding.

        Parameters
        ----------
        cov
            Covariate index.

        Returns
        -------
        Covariate name.

        """
        return "cov" + str(cov)

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

        Returns
        -------
            Tensors that can be used for inference.
            Keys:
            * ``'expr'``: Expression.
            * ``'batch'``: Batch covariate.
            * ``'cat'``: All covariates that require one-hot encoding.
                         List of tensors with dim = n_samples * 1.
                         If absent returns empty list.
            * ``'cont'``: All covariates that are already continous.
                          Includes continous and embedded
                          categorical covariates.
                          Single tensor of
                          dim = n_samples * n_concatenated_cov_features.
                          If absent returns None.
        """
        cov = self._get_cov(tensors=tensors)
        input_dict = {
            "expr": tensors[REGISTRY_KEYS.X_KEY],
            "batch": tensors[REGISTRY_KEYS.BATCH_KEY],
            "cat": cov["cat"],
            "cont": cov["cont"],
        }
        return input_dict

    def _get_inference_cycle_input(
        self,
        tensors: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        selected_batch: torch.Tensor,
        **kwargs,
    ) -> dict[str, torch.Tensor | list[torch.Tensor] | None]:
        """In. tensors, gen. outputs, and cycle batch -> cycle inference inputs.

        Parameters
        ----------
        tensors
            Input tensors.
        generative_outputs
            Outputs of the generative pass.
        selected_batch
            Batch covariate to be used for the cycle inference.
            dim = n_samples * 1
        kwargs
            Not used. Added for inheritance compatibility.

        Returns
        -------
        Tensors that can be used for cycle inference.
        Keys:
        * ``'expr'``: Expression.
        * ``'batch'``: Batch covariate.
        * ``'cat'``: All covariates that require one-hot encoding.
                     List of tensors with dim = n_samples * 1.
                     If absent returns empty list.
        * ``'cont'``: All covariates that are already continous.
                      Includes continous and embedded categorical covariates.
                      Single tensor of
                      dim = n_samples * n_concatenated_cov_features.
                      If absent returns None.

        Note: cycle covariates differ from the original publication.
        Instead of mock covariates the real input covaiates are used in cycle.

        """
        cov = self._get_cov(tensors=tensors)
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
        """In. tensors, inf. outputs, and cycle batch info -> generative inputs.

        Parameters
        ----------
        tensors
            Input tensors.
        inference_outputs
            Outputs of the inference pass.
        selected_batch
            Batch covariate to be used for the cycle expression generation.
            dim = n_samples * 1
        kwargs
            Not used. Added for inheritance compatibility.

        Returns
        -------
        Tensors that can be used for normal and cycle generation.
        Keys:
        * ``'z'``: Latent representation.
        * ``'batch'``: Batch covariates.
                       Dict with keys ``'x'`` for normal and
                       ``'y'`` for cycle pass.
        * ``'cat'``: All covariates that require one-hot encoding.
                     List of tensors with dim = n_samples * 1.
                     If absent returns empty list.
        * ``'cont'``: All covariates that are already continous.
                      Includes continous and embedded categorical covariates.
                      Single tensor of
                      dim = n_samples * n_concatenated_cov_features.
                      If absent returns None.

        Note: cycle covariates differ from the original publication.
        Instead of mock covariates the real input covaiates are used in cycle.

        """
        z = inference_outputs["z"]
        cov = self._get_cov(tensors=tensors)
        batch = {"x": tensors["batch"], "y": selected_batch}

        input_dict = {"z": z, "batch": batch, "cat": cov["cat"], "cont": cov["cont"]}
        return input_dict

    @auto_move_data  # TODO remove?
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
                      Single tensor of
                      dim = n_samples * n_concatenated_cov_features.
                      If absent returns None.

        """
        cat_parts = []
        cont_parts = []
        if REGISTRY_KEYS.CONT_COVS_KEY in tensors:
            cont_parts.append(tensors[REGISTRY_KEYS.CONT_COVS_KEY])
        if REGISTRY_KEYS.CAT_COVS_KEY in tensors:
            cat = torch.split(tensors[REGISTRY_KEYS.CAT_COVS_KEY], 1, dim=1)
            if self.embed_cat:
                for idx, tensor in enumerate(cat):
                    cont_parts.append(self.compute_embedding(self._cov_idx_name(idx), tensor))
            else:
                cat_parts.extend(cat)
        cov = {
            "cat": cat_parts,
            "cont": torch.concat(cont_parts, dim=1) if len(cont_parts) > 0 else None,
        }
        return cov

    @staticmethod
    def _merge_batch_cov(
        cat: list[torch.Tensor],
        batch: torch.Tensor,
    ) -> list[torch.Tensor]:
        """Merge batch and continuous covs for input into encoder and decoder.

        Parameters
        ----------
        cat
            Categorical covariates.
            List of tensors with dim = n_samples * 1.
        batch
            Batch covariate.
            dim = n_samples * 1

        Returns
        -------
        Single list with batch and categorical covariates.

        """
        return [batch] + cat

    @auto_move_data
    def inference(
        self,
        expr: torch.Tensor,
        batch: torch.Tensor,
        cat: list[torch.Tensor],
        cont: torch.Tensor | None,
    ) -> dict[str, torch.Tensor]:
        """Inference: expression & cov -> latent representation.

        Parameters
        ----------
        expr
            Expression data.
        batch
            Batch covariate.
        cat
            All covariates that require one-hot encoding.
        cont
            All covariates that are already continous.
            Includes continous and embedded categorical covariates.

        Returns
        -------
        Predicted mean (``'z_m'``) and variance (``'z_v'``)
        of the latent distribution as wll as a sample (``'z'``) from it.

        """
        z = self.encoder(x=expr, cat_list=self._merge_batch_cov(cat=cat, batch=batch), cont=cont)
        return {"z": z["y"], "z_m": z["y_m"], "z_v": z["y_v"]}

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        batch: dict[str, torch.Tensor],
        cat: list[torch.Tensor],
        cont: torch.Tensor | None,
        x_x: bool = True,
        x_y: bool = True,
    ) -> dict[str, torch.Tensor]:
        """Generation: latent representation & cov -> expression.

        Parameters
        ----------
        z
            Latent representation.
        batch
            Batch covariate for normal (``'x'``) and cycle (``'y'``) generation.
        cat
            All covariates that require one-hot encoding.
        cont
            All covariates that are already continous.
            Includes continous and embedded categorical covariates.
        x_x
            Decode to original batch.
        x_y
            Decode to cycle batch.

        Returns
        -------
        Predicted mean (``'x_m'``) and variance (``'x_v'``)
        of the expression distribution as wll as a sample (``'x'``) from it.
        Same outputs are returned for the cycle generation with ``'x'``
        in keys being replaced by ``'y'``.
        """

        def outputs(
            name: str,
            res: dict,
            x: torch.Tensor,
            batch: torch.Tensor,
            cat: list[torch.Tensor],
            cont: torch.Tensor | None,
        ):
            """Helper to compute generative outputs for normal and cycle pass.

            Adds generative outputs directly to the ``res`` dict.

            Parameters
            ----------
            name
                Name prepended to the keys added to the ``res`` dict.
            res
                Dict to store generative outputs in.
                Mean is stored in ``'name_m'``, variance to ``'name_v'``
                and sample to ``'name'``.
            x
                Latent representation.
            batch
                Batch covariate.
            cat
                All covariates that require one-hot encoding.
            cont
                All covariates that are already continous.
                Includes continous and embedded categorical covariates.
            """
            res_sub = self.decoder(
                x=x, cat_list=self._merge_batch_cov(cat=cat, batch=batch), cont=cont
            )
            res[name] = res_sub["y"]
            res[name + "_m"] = res_sub["y_m"]
            res[name + "_v"] = res_sub["y_v"]

        res = {}
        if x_x:
            outputs(name="x", res=res, x=z, batch=batch["x"], cat=cat, cont=cont)
        if x_y:
            outputs(name="y", res=res, x=z, batch=batch["y"], cat=cat, cont=cont)
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
        """Forward pass through the network.

        Parameters
        ----------
        tensors
            Input tensors.
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
            Whether to compute loss on forward pass.

        Returns
        -------
        Inference outputs, generative outputs of the normal pass,
        and optionally loss components.
        Inference normal and cycle outputs are combined into a single dict.
        Thus, the keys of cycle inference outputs are modified by replacing
        ``'z'`` with ``'z_cyc'``.
        """
        # TODO could disable computation of cycle if cycle loss
        #  will not be used (weight = 0).
        #  Cycle loss is not expected to be disabled in practice
        #  for typical use cases.

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

        # Combine outputs of all forward pass components
        # (first and cycle pass) into a single dict,
        # separately for inference and generative outputs
        # Rename keys in outputs of cycle pass
        # to be distinguishable from the first pass
        # for the merging into a single dict
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
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        kl_weight: float = 1.0,
        reconstruction_weight: float = 1.0,
        z_distance_cycle_weight: float = 2.0,
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
            generative_outputs["x_m"], x_true, generative_outputs["x_v"]
        ).sum(dim=1)
        reconst_loss = reconst_loss_x

        # Kl divergence on latent space
        kl_divergence_z = self.prior.kl(
            m_q=inference_outputs["z_m"],
            v_q=inference_outputs["z_v"],
            z=inference_outputs["z"],
        )

        def z_dist(
            z_x_m: torch.Tensor,
            z_y_m: torch.Tensor,
        ) -> torch.Tensor:
            """MSE loss between standardised inputs.

            MSE loss should be computed on standardized latent representations
            as else model can learn to cheat the MSE loss
            by setting the latent representations to smaller numbers.
            Standardizer is fitted on concatenation of both inputs.

            Parameters
            ----------
            z_x_m
                First input.
            z_y_m
                Second input.

            Returns
            -------
            The loss.
            dim = n_samples * 1
            """
            # Standardise data (jointly both z-s) before MSE calculation
            z = torch.concat([z_x_m, z_y_m])
            means = z.mean(dim=0, keepdim=True)
            stds = z.std(dim=0, keepdim=True)

            def standardize(x: torch.Tensor) -> torch.Tensor:
                """Helper function to standardize a tensor.

                Mean and variance from the outer scope are used for standardization.

                Parameters
                ----------
                x
                    Input tensor.

                Returns
                -------
                Standardized tensor.
                """
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
        """For each cell randomly selects new batch different from the real one.

        Parameters
        ----------
        batch
            Real batch information for each cell.

        Returns
        -------
        Newly selected batch for each cell.
        """
        # Get available batches -
        # those that are zero will become nonzero and vice versa
        batch = torch.nn.functional.one_hot(batch.squeeze(-1), self.n_batch)
        available_batches = 1 - batch
        # Get nonzero indices for each cell -
        # batches that differ from the real batch and are thus available
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

        return randomly_selected_indices

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
