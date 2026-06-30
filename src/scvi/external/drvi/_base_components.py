from __future__ import annotations

import math
from typing import TYPE_CHECKING

import torch
from stacked_linear import StackedLinearLayer
from torch import nn

from scvi.nn import FCLayers

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal


class SplitFCLayers(FCLayers):
    """Fully-connected layers that process several latent splits in parallel.

    A thin subclass of :class:`~scvi.nn.FCLayers` for DRVI's additive decoder. The input is a 3D
    tensor of shape ``(n_obs, n_split, n_features)`` (one independent "split" per channel) instead
    of the usual 2D ``(n_obs, n_features)``. Only small construction/application seams of
    :class:`~scvi.nn.FCLayers` are overridden; ``forward`` and ``set_online_update_hooks`` (the
    scArches hooks) are inherited unchanged.

    Parameters
    ----------
    n_split
        Number of parallel splits (channels) carried in the second tensor dimension.
    reuse_weights
        Per-layer weight-sharing across splits: either a single bool applied to every layer, or a
        sequence of bools with one entry per layer (length ``n_layers``). ``True`` → that layer is
        a single shared :class:`~torch.nn.Linear` (broadcast over the split dimension); ``False``
        → per-split weights via :class:`stacked_linear.StackedLinearLayer`.
    **kwargs
        Keyword arguments for :class:`~scvi.nn.FCLayers`.
    """

    def __init__(self, *args, n_split: int = 1, reuse_weights: Sequence[bool] | bool, **kwargs):
        n_layers = kwargs.get("n_layers", 1)
        if isinstance(reuse_weights, bool):
            reuse_weights = [reuse_weights] * n_layers
        elif (
            not hasattr(reuse_weights, "__len__")
            or len(reuse_weights) != n_layers
            or not all(isinstance(r, bool) for r in reuse_weights)
        ):
            raise ValueError(
                "`reuse_weights` must be a bool or a sequence of booleans with one entry per "
                f"layer (n_layers={n_layers}), got {reuse_weights!r}."
            )
        # set before super().__init__() so the overridden _build_layer can read them
        self.n_split = n_split
        self.reuse_weights = reuse_weights
        # Update some defaults
        kwargs = {
            **kwargs,
            "use_batch_norm": False,
            "use_layer_norm": True,
        }
        super().__init__(*args, **kwargs)

    def _build_linear(self, n_in: int, n_out: int, layer_num: int) -> nn.Module:
        if self.reuse_weights[layer_num]:
            # a shared nn.Linear already broadcasts over the (n_obs, n_split, .) split dimension
            return nn.Linear(n_in, n_out, bias=self.bias)
        return StackedLinearLayer(self.n_split, n_in, n_out, bias=self.bias)

    def _build_layer(self, n_in: int, n_out: int, layer_num: int) -> nn.Module:
        return nn.Sequential(
            self._build_linear(
                n_in + self.n_cov * self.inject_into_layer(layer_num),
                n_out,
                layer_num,
            ),
            nn.BatchNorm1d(self.n_split * n_out, momentum=0.01, eps=0.001)
            if self.use_batch_norm
            else None,
            # Important: layer norm is applied independently per split, so
            # LayerNorm(n_split, n_out) is not correct.
            nn.LayerNorm(n_out, elementwise_affine=False) if self.use_layer_norm else None,
            self.activation_fn() if self.use_activation else None,
            nn.Dropout(p=self.dropout_rate) if self.dropout_rate > 0 else None,
        )

    def _is_linear_layer(self, layer: nn.Module) -> bool:
        return isinstance(layer, (nn.Linear, StackedLinearLayer))

    def _apply_layer(self, layer, x, cov_list, layer_index):
        if self._is_linear_layer(layer) and self.inject_into_layer(layer_index):
            if x.dim() not in (3, 4):
                raise ValueError(
                    "SplitFCLayers expects a 3D (n_obs, n_split, n_features) input — or 4D with a "
                    f"leading n_samples axis — when injecting covariates, got a {x.dim()}D tensor."
                )
            # broadcast each covariate (n_obs, o_dim) over the split (and any leading n_samples)
            cov_list_layer = [o.unsqueeze(-2).expand(*x.shape[:-1], o.size(-1)) for o in cov_list]
            if cov_list_layer:
                x = torch.cat((x, *cov_list_layer), dim=-1)
        return layer(x)

    def _apply_batch_norm(self, layer, x):
        # batch norm over n_split * n_hidden features: fold all leading dims into the batch axis
        return layer(x.reshape(-1, x.shape[-2] * x.shape[-1])).reshape(x.shape)


class DecoderDRVI(nn.Module):
    """DRVI additive decoder: split the latent, decode each split, then aggregate.

    Reuses scvi's head structure (``px_scale_decoder``, ``px_r_decoder``, ``px_dropout_decoder``)
    but, unlike :class:`~scvi.nn.DecoderSCVI`, returns the per-gene parameters in **log space**
    (aggregated scale logits, dispersion logits and zero-inflation logits). The library-size,
    softmax and exp transforms into count space are applied by the module (see
    :meth:`~scvi.external.drvi.DRVIModule.generative`).

    The latent ``z`` of dimension ``n_latent`` is mapped into ``n_split`` independent groups, each
    decoded by :class:`SplitFCLayers`, and the per-split parameters are aggregated over the split
    dimension (in log space; the count-space transforms are applied by the module).

    Parameters
    ----------
    n_input
        Dimensionality of the latent space (``n_latent``); the dimension that is split.
    n_output
        Number of genes.
    n_split
        Number of latent splits. Must divide ``n_input`` for both ``split_mask`` and ``split_map``.
    split_method
        How the latent is mapped to splits.

        * ``"split_mask"`` — reshape into ``n_split`` contiguous chunks and place each chunk on
          its own split, zeroing the other chunks. Each split keeps the full latent width
          ``n_input``; e.g. with ``n_input=10`` and ``n_split=2`` the latent
          ``[1..10]`` becomes ``[[1,2,3,4,5,0,0,0,0,0], [0,0,0,0,0,6,7,8,9,10]]``.
        * ``"split_map"`` — reshape into ``n_split`` chunks of size ``n_input // n_split``
          and apply a learned per-split linear map (:class:`stacked_linear.StackedLinearLayer`).
    n_split_output
        Per-split projection output width, i.e. the input dimension to each split's decoder body.
        ``"auto"`` (default) uses ``n_input`` (= ``n_latent``). For ``"split_map"`` it is the
        output size of the learned per-split projection; for ``"split_mask"`` it must equal
        ``n_input``.
    split_aggregation
        How per-split parameters are combined over the split dimension.

        * ``"mean"`` — ``sum / n_split``.
        * ``"logsumexp"`` — ``logsumexp - log(n_split)`` (additive decoder in log space).
    n_continuous_cov
        Number of continuous covariates injected into the decoder layers.
    reuse_weights
        Weight-sharing policy across splits, applied to the FC body (its ``n_layers`` hidden
        layers) and the parameter heads ("last"):

        * ``"everywhere"`` (default) — hidden layers and heads all shared.
        * ``"hidden"`` — hidden layers shared, heads per-split.
        * ``"last"`` — hidden layers per-split, heads shared.
        * ``"hidden_except_first"`` — first hidden layer per-split, the rest and heads shared.
        * ``"nowhere"`` — hidden layers and heads all per-split.

        Per-split weights use :class:`stacked_linear.StackedLinearLayer`; shared weights are a
        single :class:`~torch.nn.Linear` broadcast over the split dimension.
    model_cell_dispersion
        Whether to build a per-cell dispersion head (``px_r_decoder``). If ``False`` (default), the
        decoder returns ``None`` for the dispersion and the module uses a shared dispersion
        parameter instead (needed only for ``dispersion="gene-cell"``).
    model_zero_inflation
        Whether to build a zero-inflation head (``px_dropout_decoder``). If ``False`` (default),
        the decoder returns ``None`` for the zero-inflation logits (needed only for the ``zinb``
        likelihood).
    **kwargs
        Keyword arguments for :class:`SplitFCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_split: int,
        split_method: Literal["split_mask", "split_map"] = "split_map",
        n_split_output: int | Literal["auto"] = "auto",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        n_cat_list: Sequence[int] | None = None,
        n_continuous_cov: int = 0,
        n_layers: int = 1,
        n_hidden: int = 128,
        reuse_weights: Literal[
            "everywhere", "last", "hidden", "nowhere", "hidden_except_first"
        ] = "everywhere",
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        use_layer_norm: bool = True,
        dropout_rate: float = 0.0,
        model_cell_dispersion: bool = False,
        model_zero_inflation: bool = False,
        **kwargs,
    ):
        super().__init__()
        if n_input % n_split != 0:
            raise ValueError(f"`{split_method}` requires `n_latent` divisible by `n_split`.")

        self.n_latent = n_input
        self.n_split = n_split
        self.split_method = split_method
        self.split_aggregation = split_aggregation
        self.inspect_mode = False

        if n_split_output == "auto":
            n_split_output = self.n_latent
        self.n_split_output = n_split_output

        if self.split_method == "split_mask":
            if n_split_output != self.n_latent:
                raise ValueError(
                    "For split_mask, n_split_output must be 'auto' or equal to n_latent "
                    f"({self.n_latent}), got {n_split_output}."
                )
        elif self.split_method == "split_map":
            self.split_transform = StackedLinearLayer(
                self.n_split, self.n_latent // self.n_split, n_split_output, bias=False
            )

        # Map the weight-sharing policy to per-layer decisions for the FC body (hidden layers)
        # and the parameter heads ("last"). See SplitFCLayers for the per-layer reuse mechanism.
        if reuse_weights in ("everywhere", "hidden"):
            hidden_reuse = [True] * n_layers
        elif reuse_weights == "hidden_except_first":
            hidden_reuse = [i != 0 for i in range(n_layers)]
        else:  # "last", "nowhere"
            hidden_reuse = [False] * n_layers
        last_reuse = reuse_weights in ("everywhere", "last", "hidden_except_first")

        # per-split decoder body operating on (*, n_split, n_split_output)
        self.px_decoder = SplitFCLayers(
            n_in=n_split_output,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_cont=n_continuous_cov,
            n_layers=n_layers,
            n_hidden=n_hidden,
            n_split=n_split,
            reuse_weights=hidden_reuse,
            dropout_rate=dropout_rate,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        def _make_head(n_out: int) -> nn.Module:
            if last_reuse:
                return nn.Linear(n_hidden, n_out)
            return StackedLinearLayer(n_split, n_hidden, n_out)

        self.px_scale_decoder = _make_head(n_output)
        self.px_r_decoder = _make_head(n_output) if model_cell_dispersion else None
        self.px_dropout_decoder = _make_head(n_output) if model_zero_inflation else None

    def _apply_split(self, z: torch.Tensor) -> torch.Tensor:
        """Map latent ``(*, n_latent)`` to splits ``(*, n_split, n_split_output)``."""
        *lead_shape, n_latent = z.shape
        if self.split_method == "split_mask":
            zt = z.reshape(*lead_shape, self.n_split, n_latent // self.n_split).transpose(-2, -1)
            blocks = torch.diag_embed(zt)
            return blocks.reshape(*lead_shape, self.n_split, n_latent)
        elif self.split_method == "split_map":
            z = z.reshape(*lead_shape, self.n_split, n_latent // self.n_split)
            return self.split_transform(z)
        else:
            raise ValueError(f"Invalid split_method: {self.split_method}")

    def _aggregate(self, x: torch.Tensor) -> torch.Tensor:
        """Aggregate per-split params ``(*, n_split, n_genes)`` over the split dimension."""
        n_split = x.shape[-2]
        if self.split_aggregation == "logsumexp":
            return torch.logsumexp(x, dim=-2) - math.log(n_split)
        elif self.split_aggregation == "mean":
            return x.sum(dim=-2) / n_split
        else:
            raise ValueError(f"Invalid split_aggregation: {self.split_aggregation}")

    def forward(
        self,
        z: torch.Tensor,
        *cat_list: int,
        cont: torch.Tensor | None = None,
    ):
        """Decode ``z`` into **log-space** per-gene parameters.

        Returns ``(px_scale_logit, px_r_logit, px_dropout_logit, px_scale_logit_per_split)``: the
        aggregated per-gene scale logits (log space, before softmax), the per-cell dispersion
        logits, the zero-inflation logits, and the per-split scale logits before aggregation
        (only when ``inspect_mode`` is set, else ``None``). The dispersion and zero-inflation
        outputs are ``None`` unless the corresponding head was built (``model_cell_dispersion`` /
        ``model_zero_inflation``). The library-size, softmax and exp transforms are applied by the
        module (not here). Any number of leading dimensions (e.g. an ``n_samples`` axis) is
        supported and preserved: the split transform, the per-split FC layers and the aggregation
        all act on the last one or two dimensions.
        """
        z_split = self._apply_split(z)  # (*, n_split, n_split_output)
        h = self.px_decoder(z_split, *cat_list, cont=cont)  # (*, n_split, n_hidden)

        # per-split scale logits aggregated over splits, kept in log space
        px_scale_logit_per_split = self.px_scale_decoder(h)  # (*, n_split, n_genes)
        px_scale_logit = self._aggregate(px_scale_logit_per_split)  # (*, n_genes)
        px_dropout_logit, px_r_logit = None, None
        if self.px_r_decoder is not None:
            px_r_logit = self._aggregate(self.px_r_decoder(h))  # heuristic
        if self.px_dropout_decoder is not None:
            px_dropout_logit = -self._aggregate(-self.px_dropout_decoder(h))  # heuristic

        if not self.inspect_mode:
            px_scale_logit_per_split = None

        return px_scale_logit, px_r_logit, px_dropout_logit, px_scale_logit_per_split
