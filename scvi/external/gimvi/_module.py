"""Main module."""
from typing import List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Normal, Poisson
from torch.distributions import kl_divergence as kl
from torch.nn import ModuleList

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder, MultiDecoder, MultiEncoder, one_hot

torch.backends.cudnn.benchmark = True


class JVAE(BaseModuleClass):
    """Joint variational auto-encoder for imputing missing genes in spatial data.

    Implementation of gimVI :cite:p:`Lopez19`.

    Parameters
    ----------
    dim_input_list
        List of number of input genes for each dataset. If
            the datasets have different sizes, the dataloader will loop on the
            smallest until it reaches the size of the longest one
    total_genes
        Total number of different genes
    indices_mappings
        list of mapping the model inputs to the model output
        Eg: ``[[0,2], [0,1,3,2]]`` means the first dataset has 2 genes that will be reconstructed at location ``[0,2]``
        the second dataset has 4 genes that will be reconstructed at ``[0,1,3,2]``
    gene_likelihoods
        list of distributions to use in the generative process 'zinb', 'nb', 'poisson'
    model_library_bools bool list
        model or not library size with a latent variable or use observed values
    library_log_means np.ndarray list
        List of 1 x n_batch array of means of the log library sizes.
        Parameterizes prior on library size if not using observed library sizes.
    library_log_vars np.ndarray list
        List of 1 x n_batch array of variances of the log library sizes.
        Parameterizes prior on library size if not using observed library sizes.
    n_latent
        dimension of latent space
    n_layers_encoder_individual
        number of individual layers in the encoder
    n_layers_encoder_shared
        number of shared layers in the encoder
    dim_hidden_encoder
        dimension of the hidden layers in the encoder
    n_layers_decoder_individual
        number of layers that are conditionally batchnormed in the encoder
    n_layers_decoder_shared
        number of shared layers in the decoder
    dim_hidden_decoder_individual
        dimension of the individual hidden layers in the decoder
    dim_hidden_decoder_shared
        dimension of the shared hidden layers in the decoder
    dropout_rate_encoder
        dropout encoder
    dropout_rate_decoder
        dropout decoder
    n_batch
        total number of batches
    n_labels
        total number of labels
    dispersion
        See ``vae.py``
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.

    """

    def __init__(
        self,
        dim_input_list: List[int],
        total_genes: int,
        indices_mappings: List[Union[np.ndarray, slice]],
        gene_likelihoods: List[str],
        model_library_bools: List[bool],
        library_log_means: List[Optional[np.ndarray]],
        library_log_vars: List[Optional[np.ndarray]],
        n_latent: int = 10,
        n_layers_encoder_individual: int = 1,
        n_layers_encoder_shared: int = 1,
        dim_hidden_encoder: int = 64,
        n_layers_decoder_individual: int = 0,
        n_layers_decoder_shared: int = 0,
        dim_hidden_decoder_individual: int = 64,
        dim_hidden_decoder_shared: int = 64,
        dropout_rate_encoder: float = 0.2,
        dropout_rate_decoder: float = 0.2,
        n_batch: int = 0,
        n_labels: int = 0,
        dispersion: str = "gene-batch",
        log_variational: bool = True,
    ):
        super().__init__()

        self.n_input_list = dim_input_list
        self.total_genes = total_genes
        self.indices_mappings = indices_mappings
        self.gene_likelihoods = gene_likelihoods
        self.model_library_bools = model_library_bools
        for mode in range(len(dim_input_list)):
            if self.model_library_bools[mode]:
                self.register_buffer(
                    f"library_log_means_{mode}",
                    torch.from_numpy(library_log_means[mode]).float(),
                )
                self.register_buffer(
                    f"library_log_vars_{mode}",
                    torch.from_numpy(library_log_vars[mode]).float(),
                )

        self.n_latent = n_latent

        self.n_batch = n_batch
        self.n_labels = n_labels

        self.dispersion = dispersion
        self.log_variational = log_variational

        self.z_encoder = MultiEncoder(
            n_heads=len(dim_input_list),
            n_input_list=dim_input_list,
            n_output=self.n_latent,
            n_hidden=dim_hidden_encoder,
            n_layers_individual=n_layers_encoder_individual,
            n_layers_shared=n_layers_encoder_shared,
            dropout_rate=dropout_rate_encoder,
            return_dist=True,
        )

        self.l_encoders = ModuleList(
            [
                Encoder(
                    self.n_input_list[i],
                    1,
                    n_layers=1,
                    dropout_rate=dropout_rate_encoder,
                    return_dist=True,
                )
                if self.model_library_bools[i]
                else None
                for i in range(len(self.n_input_list))
            ]
        )

        self.decoder = MultiDecoder(
            self.n_latent,
            self.total_genes,
            n_hidden_conditioned=dim_hidden_decoder_individual,
            n_hidden_shared=dim_hidden_decoder_shared,
            n_layers_conditioned=n_layers_decoder_individual,
            n_layers_shared=n_layers_decoder_shared,
            n_cat_list=[self.n_batch],
            dropout_rate=dropout_rate_decoder,
        )

        if self.dispersion == "gene":
            self.px_r = torch.nn.Parameter(torch.randn(self.total_genes))
        elif self.dispersion == "gene-batch":
            self.px_r = torch.nn.Parameter(torch.randn(self.total_genes, n_batch))
        elif self.dispersion == "gene-label":
            self.px_r = torch.nn.Parameter(torch.randn(self.total_genes, n_labels))
        else:  # gene-cell
            pass

    def sample_from_posterior_z(
        self, x: torch.Tensor, mode: int = None, deterministic: bool = False
    ) -> torch.Tensor:
        """Sample tensor of latent values from the posterior.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input)``
        mode
            head id to use in the encoder
        deterministic
            bool - whether to sample or not

        Returns
        -------
        type
            tensor of shape ``(batch_size, n_latent)``
        """
        if mode is None:
            if len(self.n_input_list) == 1:
                mode = 0
            else:
                raise Exception("Must provide a mode when having multiple datasets")
        outputs = self.inference(x, mode)
        qz_m = outputs["qz"].loc
        z = outputs["z"]
        if deterministic:
            z = qz_m
        return z

    def sample_from_posterior_l(
        self, x: torch.Tensor, mode: int = None, deterministic: bool = False
    ) -> torch.Tensor:
        """Sample the tensor of library sizes from the posterior.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input)``
            or ``(batch_size, n_input_fish)`` depending on the mode
        mode
            head id to use in the encoder
        deterministic
            bool - whether to sample or not

        Returns
        -------
        type
            tensor of shape ``(batch_size, 1)``
        """
        inference_out = self.inference(x, mode)
        return (
            inference_out["ql"].loc
            if (deterministic and inference_out["ql"] is not None)
            else inference_out["library"]
        )

    def sample_scale(
        self,
        x: torch.Tensor,
        mode: int,
        batch_index: torch.Tensor,
        y: Optional[torch.Tensor] = None,
        deterministic: bool = False,
        decode_mode: Optional[int] = None,
    ) -> torch.Tensor:
        """Return the tensor of predicted frequencies of expression.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input)``
            or ``(batch_size, n_input_fish)`` depending on the mode
        mode
            int encode mode (which input head to use in the model)
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        y
            tensor of cell-types labels with shape ``(batch_size, n_labels)``
        deterministic
            bool - whether to sample or not
        decode_mode
            int use to a decode mode different from encoding mode

        Returns
        -------
        type
            tensor of predicted expression
        """
        gen_out = self._run_forward(
            x,
            mode,
            batch_index,
            y=y,
            deterministic=deterministic,
            decode_mode=decode_mode,
        )
        return gen_out["px_scale"]

    # This is a potential wrapper for a vae like get_sample_rate
    def get_sample_rate(self, x, batch_index, *_, **__):
        """Get the sample rate for the model."""
        return self.sample_rate(x, 0, batch_index)

    def _run_forward(
        self,
        x: torch.Tensor,
        mode: int,
        batch_index: torch.Tensor,
        y: Optional[torch.Tensor] = None,
        deterministic: bool = False,
        decode_mode: int = None,
    ) -> dict:
        """Run the forward pass of the model."""
        if decode_mode is None:
            decode_mode = mode
        inference_out = self.inference(x, mode)
        if deterministic:
            z = inference_out["qz"].loc
            if inference_out["ql"] is not None:
                library = inference_out["ql"].loc
            else:
                library = inference_out["library"]
        else:
            z = inference_out["z"]
            library = inference_out["library"]
        gen_out = self.generative(z, library, batch_index, y, decode_mode)
        return gen_out

    def sample_rate(
        self,
        x: torch.Tensor,
        mode: int,
        batch_index: torch.Tensor,
        y: Optional[torch.Tensor] = None,
        deterministic: bool = False,
        decode_mode: int = None,
    ) -> torch.Tensor:
        """Returns the tensor of scaled frequencies of expression.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input)``
            or ``(batch_size, n_input_fish)`` depending on the mode
        y
            tensor of cell-types labels with shape ``(batch_size, n_labels)``
        mode
            int encode mode (which input head to use in the model)
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        deterministic
            bool - whether to sample or not
        decode_mode
            int use to a decode mode different from encoding mode

        Returns
        -------
        type
            tensor of means of the scaled frequencies
        """
        gen_out = self._run_forward(
            x,
            mode,
            batch_index,
            y=y,
            deterministic=deterministic,
            decode_mode=decode_mode,
        )
        return gen_out["px_rate"]

    def reconstruction_loss(
        self,
        x: torch.Tensor,
        px_rate: torch.Tensor,
        px_r: torch.Tensor,
        px_dropout: torch.Tensor,
        mode: int,
    ) -> torch.Tensor:
        """Compute the reconstruction loss."""
        reconstruction_loss = None
        if self.gene_likelihoods[mode] == "zinb":
            reconstruction_loss = (
                -ZeroInflatedNegativeBinomial(
                    mu=px_rate, theta=px_r, zi_logits=px_dropout
                )
                .log_prob(x)
                .sum(dim=-1)
            )
        elif self.gene_likelihoods[mode] == "nb":
            reconstruction_loss = (
                -NegativeBinomial(mu=px_rate, theta=px_r).log_prob(x).sum(dim=-1)
            )
        elif self.gene_likelihoods[mode] == "poisson":
            reconstruction_loss = -Poisson(px_rate).log_prob(x).sum(dim=1)
        return reconstruction_loss

    def _get_inference_input(self, tensors):
        """Get the input for the inference model."""
        return {"x": tensors[REGISTRY_KEYS.X_KEY]}

    def _get_generative_input(self, tensors, inference_outputs):
        """Get the input for the generative model."""
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        return {"z": z, "library": library, "batch_index": batch_index, "y": y}

    @auto_move_data
    def inference(self, x: torch.Tensor, mode: Optional[int] = None) -> dict:
        """Run the inference model."""
        x_ = x
        if self.log_variational:
            x_ = torch.log(1 + x_)

        qz, z = self.z_encoder(x_, mode)
        ql, library = None, None
        if self.model_library_bools[mode]:
            ql, library = self.l_encoders[mode](x_)
        else:
            library = torch.log(torch.sum(x, dim=1)).view(-1, 1)

        return {"qz": qz, "z": z, "ql": ql, "library": library}

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        y: Optional[torch.Tensor] = None,
        mode: Optional[int] = None,
    ) -> dict:
        """Run the generative model."""
        px_scale, px_r, px_rate, px_dropout = self.decoder(
            z, mode, library, self.dispersion, batch_index, y
        )
        if self.dispersion == "gene-label":
            px_r = F.linear(one_hot(y, self.n_labels), self.px_r)
        elif self.dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r.view(1, self.px_r.size(0))
        px_r = torch.exp(px_r)

        px_scale = px_scale / torch.sum(
            px_scale[:, self.indices_mappings[mode]], dim=1
        ).view(-1, 1)
        px_rate = px_scale * torch.exp(library)

        return {
            "px_scale": px_scale,
            "px_r": px_r,
            "px_rate": px_rate,
            "px_dropout": px_dropout,
        }

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        mode: Optional[int] = None,
        kl_weight=1.0,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Return the reconstruction loss and the Kullback divergences.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input)``
            or ``(batch_size, n_input_fish)`` depending on the mode
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        y
            tensor of cell-types labels with shape (batch_size, n_labels)
        mode
            indicates which head/tail to use in the joint network


        Returns
        -------
        the reconstruction loss and the Kullback divergences
        """
        if mode is None:
            if len(self.n_input_list) == 1:
                mode = 0
            else:
                raise Exception("Must provide a mode")
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        qz = inference_outputs["qz"]
        ql = inference_outputs["ql"]
        px_rate = generative_outputs["px_rate"]
        px_r = generative_outputs["px_r"]
        px_dropout = generative_outputs["px_dropout"]

        # mask loss to observed genes
        mapping_indices = self.indices_mappings[mode]
        reconstruction_loss = self.reconstruction_loss(
            x,
            px_rate[:, mapping_indices],
            px_r[:, mapping_indices],
            px_dropout[:, mapping_indices],
            mode,
        )

        # KL Divergence
        mean = torch.zeros_like(qz.loc)
        scale = torch.ones_like(qz.scale)
        kl_divergence_z = kl(qz, Normal(mean, scale)).sum(dim=1)

        if self.model_library_bools[mode]:
            library_log_means = getattr(self, f"library_log_means_{mode}")
            library_log_vars = getattr(self, f"library_log_vars_{mode}")

            local_library_log_means = F.linear(
                one_hot(batch_index, self.n_batch), library_log_means
            )
            local_library_log_vars = F.linear(
                one_hot(batch_index, self.n_batch), library_log_vars
            )
            kl_divergence_l = kl(
                ql,
                Normal(local_library_log_means, local_library_log_vars.sqrt()),
            ).sum(dim=1)
        else:
            kl_divergence_l = torch.zeros_like(kl_divergence_z)

        kl_local = kl_divergence_l + kl_divergence_z

        loss = torch.mean(reconstruction_loss + kl_weight * kl_local) * x.size(0)

        return LossOutput(
            loss=loss, reconstruction_loss=reconstruction_loss, kl_local=kl_local
        )
