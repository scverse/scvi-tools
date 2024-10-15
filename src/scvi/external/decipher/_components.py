from collections.abc import Sequence

import numpy as np
import torch
import torch.nn as nn


class ConditionalDenseNN(nn.Module):
    """Dense neural network with multiple outputs, optionally conditioned on a context variable.

    (Derived from pyro.nn.dense_nn.ConditionalDenseNN with some modifications [1])

    Parameters
    ----------
    input_dim : int
        Dimension of the input
    hidden_dims : sequence of ints
        Dimensions of the hidden layers (excluding the output layer)
    output_dims : sequence of ints (optional)
        Dimensions of each output layer
        Default: (1,)
    context_dim : int (optional)
        Dimension of the context input.
        Default: 0. No context input.
    deep_context_injection : bool (optional)
        If True, inject the context into every hidden layer.
        If False, only inject the context into the first hidden layer
        (concatenated with the input).
        Default: False.
    activation : torch.nn.Module (optional)
        Activation function to use between hidden layers (not applied to the outputs).
        Default: torch.nn.ReLU()
    """

    def __init__(
        self,
        input_dim: int,
        hidden_dims: Sequence[int],
        output_dims: Sequence = (1,),
        context_dim: int = 0,
        deep_context_injection: bool = False,
        activation=torch.nn.ReLU(),
    ):
        super().__init__()

        self.input_dim = input_dim
        self.context_dim = context_dim
        self.hidden_dims = hidden_dims
        self.output_dims = output_dims
        self.deep_context_injection = deep_context_injection
        self.n_output_layers = len(self.output_dims)
        self.output_total_dim = sum(self.output_dims)

        # The multiple outputs are computed as a single output layer, and then split
        indices = np.concatenate(([0], np.cumsum(self.output_dims)))
        self.output_slices = [slice(s, e) for s, e in zip(indices[:-1], indices[1:], strict=False)]

        # Create masked layers
        deep_context_dim = self.context_dim if self.deep_context_injection else 0
        layers = []
        batch_norms = []
        if len(hidden_dims):
            layers.append(torch.nn.Linear(input_dim + context_dim, hidden_dims[0]))
            batch_norms.append(nn.BatchNorm1d(hidden_dims[0]))
            for i in range(1, len(hidden_dims)):
                layers.append(
                    torch.nn.Linear(hidden_dims[i - 1] + deep_context_dim, hidden_dims[i])
                )
                batch_norms.append(nn.BatchNorm1d(hidden_dims[i]))

            layers.append(
                torch.nn.Linear(hidden_dims[-1] + deep_context_dim, self.output_total_dim)
            )
        else:
            layers.append(torch.nn.Linear(input_dim + context_dim, self.output_total_dim))

        self.layers = torch.nn.ModuleList(layers)

        self.f = activation
        self.batch_norms = torch.nn.ModuleList(batch_norms)

    def forward(self, x, context=None):
        if context is not None:
            # We must be able to broadcast the size of the context over the input
            context = context.expand(x.size()[:-1] + (context.size(-1),))

        h = x
        for i, layer in enumerate(self.layers):
            if self.context_dim > 0 and (self.deep_context_injection or i == 0):
                h = torch.cat([context, h], dim=-1)
            h = layer(h)
            if i < len(self.layers) - 1:
                h = self.batch_norms[i](h)
                h = self.f(h)

        if self.n_output_layers == 1:
            return h
        else:
            h = h.reshape(list(x.size()[:-1]) + [self.output_total_dim])

            if self.n_output_layers == 1:
                return h

            else:
                return tuple([h[..., s] for s in self.output_slices])
