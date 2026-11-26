from __future__ import annotations

import collections
from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.distributions import Normal

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Literal


def _identity(x):
    return x


class ConditionalBatchNorm1d(nn.Module):
    def __init__(self, num_features: int, num_classes: int, momentum: float, eps: float):
        super().__init__()
        self.num_features = num_features
        self.bn = nn.BatchNorm1d(num_features, momentum=momentum, eps=eps, affine=False)
        self.embed = nn.Embedding(num_classes, num_features * 2)
        self.embed.weight.data[:, :num_features].normal_(1, 0.02)
        self.embed.weight.data[:, num_features:].zero_()

    def forward(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        out = self.bn(x)
        gamma, beta = self.embed(y.long().ravel()).chunk(2, 1)
        return gamma.view(-1, self.num_features) * out + beta.view(-1, self.num_features)


class ConditionalLayerNorm(nn.Module):
    def __init__(self, num_features: int, num_classes: int):
        super().__init__()
        self.num_features = num_features
        self.ln = nn.LayerNorm(num_features, elementwise_affine=False)
        self.embed = nn.Embedding(num_classes, num_features * 2)
        self.embed.weight.data[:, :num_features].normal_(1, 0.02)
        self.embed.weight.data[:, num_features:].zero_()

    def forward(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        out = self.ln(x)
        gamma, beta = self.embed(y.long().ravel()).chunk(2, 1)
        return gamma.view(-1, self.num_features) * out + beta.view(-1, self.num_features)


class FCLayersX(nn.Module):
    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_continuous: int = 0,
        n_cat_list: Iterable[int] | None = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        use_activation: bool = True,
        bias: bool = True,
        inject_covariates: bool = True,
        activation_fn: Callable[[], nn.Module] = nn.ReLU,
        conditional_norm: bool = False,
        conditional_category: int = 0,
    ):
        super().__init__()
        self.inject_covariates = inject_covariates
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]

        if n_cat_list is not None:
            self.n_cat_list = [n_cat if n_cat > 1 else 0 for n_cat in n_cat_list]
        else:
            self.n_cat_list = []
        self.n_continuous = n_continuous
        self.cond_cat = conditional_category

        if conditional_norm:
            invalid_conditional_category = (
                conditional_category >= len(self.n_cat_list)
                or self.n_cat_list[conditional_category] == 0
            )
            if invalid_conditional_category:
                raise ValueError(
                    "Conditional normalization requires a categorical covariate with more than "
                    "one category."
                )

        self.n_cov = n_continuous + sum(self.n_cat_list)
        self.fc_layers = nn.Sequential(
            collections.OrderedDict(
                [
                    (
                        f"Layer {i}",
                        nn.Sequential(
                            nn.Linear(
                                n_in + self.n_cov * self.inject_into_layer(i),
                                n_out,
                                bias=bias,
                            ),
                            ConditionalBatchNorm1d(
                                n_out,
                                self.n_cat_list[self.cond_cat],
                                momentum=0.01,
                                eps=0.001,
                            )
                            if conditional_norm and use_batch_norm
                            else nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001)
                            if use_batch_norm
                            else None,
                            ConditionalLayerNorm(n_out, self.n_cat_list[self.cond_cat])
                            if conditional_norm and use_layer_norm
                            else nn.LayerNorm(n_out, elementwise_affine=False)
                            if use_layer_norm
                            else None,
                            activation_fn() if use_activation else None,
                            nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None,
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(
                        zip(layers_dim[:-1], layers_dim[1:], strict=True)
                    )
                ]
            )
        )

    def inject_into_layer(self, layer_num: int) -> bool:
        return layer_num == 0 or (layer_num > 0 and self.inject_covariates)

    def set_online_update_hooks(self, hook_first_layer: bool = True):
        self.hooks = []

        def _hook_fn_weight(grad):
            categorical_dims = sum(self.n_cat_list)
            new_grad = torch.zeros_like(grad)
            if categorical_dims > 0:
                new_grad[:, -categorical_dims:] = grad[:, -categorical_dims:]
            return new_grad

        def _hook_fn_zero_out(grad):
            return grad * 0

        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if i == 0 and not hook_first_layer:
                    continue
                if isinstance(layer, nn.Linear):
                    if self.inject_into_layer(i):
                        w = layer.weight.register_hook(_hook_fn_weight)
                    else:
                        w = layer.weight.register_hook(_hook_fn_zero_out)
                    self.hooks.append(w)
                    b = layer.bias.register_hook(_hook_fn_zero_out)
                    self.hooks.append(b)

    def forward(
        self,
        x: torch.Tensor,
        *cat_list: int,
        cont_input: torch.Tensor | None = None,
    ) -> torch.Tensor:
        one_hot_cat_list = []
        cont_list = [cont_input] if cont_input is not None else []
        cat_list = cat_list or []

        if len(self.n_cat_list) > len(cat_list):
            raise ValueError("nb. categorical args provided doesn't match init. params.")
        if (
            self.n_continuous > 0
            and cont_input is not None
            and cont_input.shape[-1] != self.n_continuous
        ):
            raise ValueError("continuous dims provided doesn't match init. params.")

        for n_cat, cat in zip(self.n_cat_list, cat_list, strict=False):
            if n_cat and cat is None:
                raise ValueError("cat not provided while n_cat != 0 in init. params.")
            if n_cat > 1:
                if cat.size(1) != n_cat:
                    one_hot_cat = nn.functional.one_hot(cat.squeeze(-1), n_cat)
                else:
                    one_hot_cat = cat
                one_hot_cat_list.append(one_hot_cat)
        cov_list = cont_list + one_hot_cat_list

        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if layer is None:
                    continue
                if isinstance(layer, (ConditionalBatchNorm1d, ConditionalLayerNorm)):
                    if x.dim() == 3:
                        x = torch.cat(
                            [
                                layer(slice_x, cat_list[self.cond_cat]).unsqueeze(0)
                                for slice_x in x
                            ],
                            dim=0,
                        )
                    else:
                        x = layer(x, cat_list[self.cond_cat])
                elif isinstance(layer, nn.BatchNorm1d):
                    if x.dim() == 3:
                        if x.device.type == "mps":
                            x = torch.cat(
                                [(layer(slice_x.clone())).unsqueeze(0) for slice_x in x], dim=0
                            )
                        else:
                            x = torch.cat([layer(slice_x).unsqueeze(0) for slice_x in x], dim=0)
                    else:
                        x = layer(x)
                else:
                    if isinstance(layer, nn.Linear) and self.inject_into_layer(i):
                        if x.dim() == 3:
                            cov_list_layer = [
                                o.unsqueeze(0).expand((x.size(0), o.size(0), o.size(1)))
                                for o in cov_list
                            ]
                        else:
                            cov_list_layer = cov_list
                        x = torch.cat((x, *cov_list_layer), dim=-1)
                    x = layer(x)
        return x


class EncoderX(nn.Module):
    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_continuous: int = 0,
        n_cat_list: Iterable[int] | None = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        distribution: str = "normal",
        var_eps: float = 1e-4,
        var_activation: Callable | None = None,
        return_dist: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.distribution = distribution
        self.var_eps = var_eps
        self.encoder = FCLayersX(
            n_in=n_input,
            n_out=n_hidden,
            n_continuous=n_continuous,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            **kwargs,
        )
        self.mean_encoder = nn.Linear(n_hidden, n_output)
        self.var_encoder = nn.Linear(n_hidden, n_output)
        self.return_dist = return_dist

        self.z_transformation = nn.Softmax(dim=-1) if distribution == "ln" else _identity
        self.var_activation = torch.exp if var_activation is None else var_activation

    def forward(
        self,
        x: torch.Tensor,
        *cat_list: int,
        cont: torch.Tensor | None = None,
    ):
        q = self.encoder(x, *cat_list, cont_input=cont)
        q_m = self.mean_encoder(q)
        q_v = self.var_activation(self.var_encoder(q)) + self.var_eps
        dist = Normal(q_m, q_v.sqrt())
        latent = self.z_transformation(dist.rsample())
        if self.return_dist:
            return dist, latent
        return q_m, q_v, latent


class DecoderSCVIX(nn.Module):
    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_continuous: int = 0,
        n_cat_list: Iterable[int] | None = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        n_conditions_output: int = 0,
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        use_layer_norm: bool = False,
        scale_activation: Literal["softmax", "softplus"] = "softmax",
        **kwargs,
    ):
        super().__init__()
        self.n_conditions_output = n_conditions_output
        self.px_decoder = FCLayersX(
            n_in=n_input,
            n_out=n_hidden,
            n_continuous=n_continuous,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        if scale_activation == "softmax":
            px_scale_activation = nn.Softmax(dim=-1)
        elif scale_activation == "softplus":
            px_scale_activation = nn.Softplus()
        else:
            raise ValueError("`scale_activation` must be 'softmax' or 'softplus'.")

        self.px_scale_decoder = nn.Sequential(
            nn.Linear(n_hidden + n_conditions_output, n_output),
            px_scale_activation,
        )
        self.px_r_decoder = nn.Linear(n_hidden + n_conditions_output, n_output)
        self.px_dropout_decoder = nn.Linear(n_hidden + n_conditions_output, n_output)

    def forward(
        self,
        dispersion: str,
        z: torch.Tensor,
        library: torch.Tensor,
        *cat_list: int,
        cont: torch.Tensor | None = None,
        output_condition: torch.Tensor | None = None,
    ):
        px = self.px_decoder(z, *cat_list, cont_input=cont)
        if output_condition is not None and self.n_conditions_output:
            one_hot_cat = nn.functional.one_hot(
                output_condition.squeeze(-1), self.n_conditions_output
            ).float()
        else:
            one_hot_cat = torch.zeros(px.size(-2), self.n_conditions_output, device=px.device)
        if px.dim() == 3:
            one_hot_cat = one_hot_cat.unsqueeze(0).expand(px.size(0), -1, -1)
        px_cat = torch.cat([px, one_hot_cat], dim=-1)

        px_scale = self.px_scale_decoder(px_cat)
        px_dropout = self.px_dropout_decoder(px_cat)
        px_rate = torch.exp(library) * px_scale
        px_r = self.px_r_decoder(px_cat) if dispersion == "gene-cell" else None
        return px_scale, px_r, px_rate, px_dropout
