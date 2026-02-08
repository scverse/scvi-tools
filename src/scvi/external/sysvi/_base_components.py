from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

import numpy as np
import torch
from torch.nn import (
    Linear,
    Module,
    Parameter,
)

from scvi import settings
from scvi.distributions import Normal
from scvi.nn import FCLayers


class EncoderDecoder(Module):
    """Module that can be used as probabilistic encoder or decoder.

    Based on inputs and optional covariates predicts output mean and variance.

    Parameters
    ----------
    n_input
        The dimensionality of the main input.
    n_output
        The dimensionality of the output.
    n_cat_list
        A list containing the number of categories for each covariate.
    n_cont
        The dimensionality of the continuous covariates.
    n_hidden
        The number of nodes per hidden layer.
    n_layers
        The number of hidden layers.
    var_mode
        How to compute variance from model outputs,
        see :class:`~scvi.external.sysvi._base_components.VarEncoder`.
        One of the following:
        * ```'sample_feature'``` - learn variance per sample and feature.
        * ```'feature'``` - learn variance per feature, constant across samples.
    var_activation
        Function used to ensure positivity of the variance.
        Defaults to :meth:`torch.exp`.
    sample
        Return samples from predicted distribution.
    kwargs
         Passed to :class:`~scvi.external.sysvi.Layers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: list[int],
        n_cont: int,
        n_hidden: int = 256,
        n_layers: int = 3,
        var_mode: Literal["sample_feature", "feature"] = "feature",
        var_activation: Callable | None = None,
        sample: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.sample = sample

        self.decoder_y = FCLayers(
            n_in=n_input,
            n_cat_list=n_cat_list,
            n_cont=n_cont,
            n_out=n_hidden,
            n_hidden=n_hidden,
            n_layers=n_layers,
            **kwargs,
        )

        self.mean_encoder = Linear(n_hidden, n_output)
        self.var_encoder = VarEncoder(n_hidden, n_output, mode=var_mode, activation=var_activation)

    def forward(
        self,
        x: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        cont: torch.Tensor | None = None,
        cat_list: list[torch.Tensor] | None = None,
    ) -> dict[str, torch.Tensor]:
        """Forward pass.

        Parameters
        ----------
        x
            Main input (i.e. expression for encoder or
            latent embedding for decoder.).
            dim = n_samples * n_input
        batch_index
            Batch index covariate.
            dim = n_samples * 1
        cont
            Optional continuous covariates.
            dim = n_samples * n_cont
        cat_list
            List of optional categorical covariates.
            Will be one hot encoded in `~scvi.nn.FCLayers`.
            Each list entry is of dim = n_samples * 1

        Returns
        -------
        Predicted mean (``'q_m'``) and variance (``'q_v'``) and
        optionally samples (``'q'``) from normal distribution
        parametrized with the predicted parameters.
        """
        cat_list = [batch_index] + cat_list
        q_ = self.decoder_y(x, *cat_list, cont=cont)
        q_m = self.mean_encoder(q_)
        if q_m.isnan().any() or q_m.isinf().any():
            warnings.warn(
                "Predicted mean contains nan or infinity values. " + "Setting to numerical.",
                stacklevel=settings.warnings_stacklevel,
            )
            q_m = torch.nan_to_num(q_m)
        q_v = self.var_encoder(q_)

        outputs = {"q_dist": Normal(q_m, q_v.sqrt() + 1e-12)}

        if self.sample:
            outputs["q"] = outputs["q_dist"].rsample()

        return outputs


class VarEncoder(Module):
    """Encode variance (strictly positive).

    Parameters
    ----------
    n_input
        Number of input dimensions.
        Used if mode is ``'sample_feature'`` to construct a network predicting
        variance from input features.
    n_output
        Number of variances to predict, matching the desired number of features
        (e.g. latent dimensions for variational encoding or output features
        for variational decoding).
    mode
        How to compute variance.
        One of the following:
        * ``'sample_feature'`` - learn variance per sample and feature.
        * ``'feature'`` - learn variance per feature, constant across samples.
    activation
        Activation function. If empty it is set to softplus.
    exp_clip
        Perform clipping before activation to prevent inf values in e**x.
        This should be useful for any activation function using e**x,
        such as exp, softplus, etc.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        mode: Literal["sample_feature", "feature", "linear"],
        activation: Callable | None = None,
        exp_clip: bool = True,
    ):
        super().__init__()

        self.clip_exp_thr = np.log(torch.finfo(torch.get_default_dtype()).max) - 1e-4
        self.exp_clip = exp_clip
        self.mode = mode
        if self.mode == "sample_feature":
            self.encoder = Linear(n_input, n_output)
        elif self.mode == "feature":
            self.var_param = Parameter(torch.zeros(1, n_output))
        else:
            raise ValueError("Mode not recognised.")
        self.activation = torch.nn.Softplus() if activation is None else activation

    def forward(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:
        """Forward pass through model.

        Parameters
        ----------
        x
            Used to encode variance if mode is ``'sample_feature'``.
            dim = n_samples x n_input

        Returns
        -------
        Predicted variance
        dim = n_samples * 1
        """
        if self.mode == "sample_feature":
            v = self.encoder(x)
        elif self.mode == "feature":
            v = self.var_param.expand(x.shape[0], -1)  # Broadcast to input size

        # Ensure that var is strictly positive via exp -
        # Bring back to non-log scale
        # Clip to range that will not be inf after exp
        # This should be useful for any activation that uses e**x
        # such as exp, softplus, etc.
        if self.exp_clip:
            v = torch.clip(v, min=-self.clip_exp_thr, max=self.clip_exp_thr)
        v = self.activation(v)
        if v.isnan().any():
            warnings.warn(
                "Predicted variance contains nan values. Setting to 0.",
                stacklevel=settings.warnings_stacklevel,
            )
            v = torch.nan_to_num(v)

        return v
