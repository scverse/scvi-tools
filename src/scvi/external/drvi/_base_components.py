from __future__ import annotations

import math
from typing import TYPE_CHECKING

import torch
from stacked_linear import StackedLinearLayer
from torch import nn

from scvi.nn import DecoderSCVI, FCLayers

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal


class SplitFCLayers(FCLayers):
    """Fully-connected layers that process several latent splits in parallel.

    A thin subclass of :class:`~scvi.nn.FCLayers` for DRVI's additive decoder. The input is a 3D
    tensor of shape ``(n_obs, n_split, n_features)`` (one independent "split" per channel) instead
    of the usual 2D ``(n_obs, n_features)``. Only the layer-construction / per-layer application
    seams of :class:`~scvi.nn.FCLayers` are overridden — ``__init__`` and ``forward`` are inherited
    unchanged.

    Parameters
    ----------
    n_split
        Number of parallel splits (channels) carried in the second tensor dimension.
    reuse_weights
        If ``True`` (default), all splits share a single :class:`~torch.nn.Linear` per layer
        (its weights broadcast over the split dimension). If ``False``, each split gets its own
        weights via :class:`stacked_linear.StackedLinearLayer`.
    **kwargs
        Keyword arguments for :class:`~scvi.nn.FCLayers`.
    """

    def __init__(self, *args, n_split: int = 1, reuse_weights: bool = True, **kwargs):
        # set before super().__init__() so the overridden _build_linear can read them
        self._n_split = n_split
        self._reuse_weights = reuse_weights
        super().__init__(*args, **kwargs)

    def _build_linear(self, n_in: int, n_out: int, bias: bool) -> nn.Module:
        if self._n_split <= 1 or self._reuse_weights:
            # a shared nn.Linear already broadcasts over the (n_obs, n_split, .) split dimension
            return nn.Linear(n_in, n_out, bias=bias)
        return StackedLinearLayer(self._n_split, n_in, n_out, bias=bias)

    def _apply_layer(self, layer, x, cov_list, layer_index):
        # recognize the stacked linear in addition to nn.Linear; broadcast covariates across the
        # split dimension of a 3D (n_obs, n_split, n_features) tensor.
        is_linear = isinstance(layer, (nn.Linear, StackedLinearLayer))
        if is_linear and self.inject_into_layer(layer_index):
            if x.dim() == 3:
                cov_list_layer = [
                    o.unsqueeze(1).expand(o.size(0), x.size(1), o.size(-1)) for o in cov_list
                ]
            else:
                cov_list_layer = cov_list
            if cov_list_layer:
                x = torch.cat((x, *cov_list_layer), dim=-1)
        return layer(x)

    def _apply_batch_norm(self, layer, x):
        # batch-norm over (n_obs * n_split) per feature; DRVI/scvi decoders default this off and
        # use layer norm (which already works on the last dim of a 3D tensor) instead.
        if x.dim() == 3:
            n_obs, n_split, n_features = x.shape
            x = layer(x.reshape(n_obs * n_split, n_features))
            return x.reshape(n_obs, n_split, n_features)
        return layer(x)


class SplitDecoder(DecoderSCVI):
    """DRVI additive decoder: split the latent, decode each split, then aggregate.

    Subclasses :class:`~scvi.nn.DecoderSCVI`, reusing its parameter heads (``px_scale_decoder``,
    ``px_r_decoder``, ``px_dropout_decoder``) and producing the same ``(px_scale, px_r, px_rate,
    px_dropout)`` tuple, so it is a drop-in for the scvi generative process and likelihoods.

    The latent ``z`` of dimension ``n_latent`` is mapped into ``n_split`` independent groups, each
    decoded by :class:`SplitFCLayers`, and the per-split parameters are aggregated over the split
    dimension before the final activation.

    Parameters
    ----------
    n_input
        Dimensionality of the latent space (``n_latent``); the dimension that is split.
    n_output
        Number of genes.
    n_split
        Number of latent splits. With ``split_diag`` it must equal ``n_input``; with ``split_map``
        it must divide ``n_input``.
    split_method
        How the latent is mapped to splits.

        * ``"split_diag"`` — ``torch.diag_embed(z)`` so split ``i`` sees only latent dim ``i``.
        * ``"split_map"`` — reshape into ``n_split`` chunks of size ``n_input // n_split``
          and apply a learned per-split linear map (:class:`stacked_linear.StackedLinearLayer`).
    split_aggregation
        How per-split parameters are combined over the split dimension.

        * ``"mean"`` — ``sum / n_split``.
        * ``"logsumexp"`` — ``logsumexp - log(n_split)`` (additive decoder in log space).
    n_continuous_cov
        Number of continuous covariates injected into the decoder layers.
    reuse_weights
        Passed to :class:`SplitFCLayers`: share decoder weights across splits (``True``) or learn
        per-split weights (``False``).
    **kwargs
        Keyword arguments for :class:`~scvi.nn.DecoderSCVI` / :class:`SplitFCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_split: int,
        split_method: Literal["split_diag", "split_map"] = "split_map",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        n_cat_list: Iterable[int] | None = None,
        n_continuous_cov: int = 0,
        n_layers: int = 1,
        n_hidden: int = 128,
        reuse_weights: bool = True,
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
        scale_activation: Literal["softmax", "softplus"] = "softmax",
        dropout_rate: float = 0.0,
        **kwargs,
    ):
        if split_method not in ("split_diag", "split_map"):
            raise ValueError("`split_method` must be one of 'split_diag', 'split_map'.")
        if split_aggregation not in ("mean", "logsumexp"):
            raise ValueError("`split_aggregation` must be one of 'mean', 'logsumexp'.")

        self.n_latent = n_input
        self.n_split = n_split
        self.split_method = split_method
        self.split_aggregation = split_aggregation
        self.inspect_mode = False
        # holds per-split pre-aggregation scale logits when inspect_mode is on
        self._split_scale_cache: torch.Tensor | None = None

        split_in = self._split_input_dim(n_input)

        # builds the parameter heads (n_hidden -> n_output) and a throwaway px_decoder we replace.
        super().__init__(
            n_input=split_in,
            n_output=n_output,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            scale_activation=scale_activation,
            **kwargs,
        )
        # per-split decoder body operating on (n_obs, n_split, split_in)
        self.px_decoder = SplitFCLayers(
            n_in=split_in,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_cont=n_continuous_cov,
            n_layers=n_layers,
            n_hidden=n_hidden,
            n_split=n_split,
            reuse_weights=reuse_weights,
            dropout_rate=dropout_rate,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )
        if self.split_method == "split_map":
            # learned per-split projection (submodule must be created after nn.Module.__init__):
            # (n_obs, n_split, split_in) -> (n_obs, n_split, split_in)
            self.split_transform = StackedLinearLayer(self.n_split, split_in, split_in, bias=False)

    def _split_input_dim(self, n_input: int) -> int:
        """Validate the split config and return the per-split input dimension."""
        if self.split_method == "split_diag":
            if n_input != self.n_split:
                raise ValueError("`split_diag` requires `n_split == n_latent`.")
            return n_input
        # split_map
        if n_input % self.n_split != 0:
            raise ValueError("`split_map` requires `n_latent` divisible by `n_split`.")
        return n_input // self.n_split

    def _apply_split(self, z: torch.Tensor) -> torch.Tensor:
        """Map latent ``(n_obs, n_latent)`` to splits ``(n_obs, n_split, split_in)``."""
        if self.split_method == "split_diag":
            return torch.diag_embed(z)
        # split_map
        z = z.reshape(z.shape[0], self.n_split, -1)
        return self.split_transform(z)

    def _aggregate(self, x: torch.Tensor) -> torch.Tensor:
        """Aggregate per-split params ``(n_obs, n_split, n_genes)`` over the split dimension."""
        n_split = x.shape[-2]
        if self.split_aggregation == "mean":
            return x.sum(dim=-2) / n_split
        # logsumexp: cancels the n_split factor of an additive (log-space) decoder
        return torch.logsumexp(x, dim=-2) - math.log(n_split)

    def forward(
        self,
        dispersion: str,
        z: torch.Tensor,
        library: torch.Tensor,
        *cat_list: int,
        cont: torch.Tensor | None = None,
    ):
        """Decode ``z`` into ZINB/NB parameters via per-split decoding and aggregation."""
        # flatten any leading dims (e.g. n_samples, n_obs) into a single observation axis
        leading = z.shape[:-1]
        z2 = z.reshape(-1, z.shape[-1])
        n_obs = z2.shape[0]
        library2 = library.reshape(-1, library.shape[-1])

        def _match(t):
            # per-cell covariates may need tiling to match the flattened (n_samples * n_obs) axis
            if t is None:
                return None
            t = t.reshape(-1, t.shape[-1]) if t.dim() > 2 else t
            if t.shape[0] != n_obs and n_obs % t.shape[0] == 0:
                t = t.repeat(n_obs // t.shape[0], 1)
            return t

        cont2 = _match(cont)
        cat2 = [_match(c) for c in cat_list]

        z_split = self._apply_split(z2)  # (n_obs, n_split, split_in)
        h = self.px_decoder(z_split, *cat2, cont=cont2)  # (n_obs, n_split, n_hidden)

        # per-split scale logits, aggregate over splits (raw, log-space for logsumexp), then the
        # aggregated parameter is either passed through the scale activation (scvi likelihoods) or
        # consumed directly in log space by the module (e.g. the parametrized/log NB and normal).
        scale_logits = self.px_scale_decoder[0](h)  # (n_obs, n_split, n_genes)
        px_agg = self._aggregate(scale_logits)
        px_scale = self.px_scale_decoder[1](px_agg)
        px_rate = torch.exp(library2) * px_scale
        # dropout / cell-wise dispersion are auxiliary: average their per-split logits
        px_dropout = self.px_dropout_decoder(h).mean(dim=-2)
        px_r = self.px_r_decoder(h).mean(dim=-2) if dispersion == "gene-cell" else None

        if self.inspect_mode:
            self._split_scale_cache = scale_logits.reshape(*leading, self.n_split, -1)

        def _unflatten(t):
            return None if t is None else t.reshape(*leading, t.shape[-1])

        return (
            _unflatten(px_scale),
            _unflatten(px_r),
            _unflatten(px_rate),
            _unflatten(px_dropout),
            _unflatten(px_agg),
        )
