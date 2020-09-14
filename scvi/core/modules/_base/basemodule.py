from abc import abstractmethod
from typing import Tuple

import torch
import torch.nn as nn


class AbstractVAE(nn.Module):
    def __init__(
        self,
    ):
        super().__init__()

    def forward(
        self,
        tensors,
        inference_kwargs={},
        generative_kwargs={},
        loss_kwargs={},
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass through the network

        Parameters
        ----------
        tensors
            tensors to pass through
        """
        inference_kwargs = dict(inference_kwargs)
        generative_kwargs = dict(generative_kwargs)
        loss_kwargs = dict(loss_kwargs)

        inference_outputs = self.inference(tensors, **inference_kwargs)
        generative_outputs = self.generative(
            tensors, inference_outputs, **generative_kwargs
        )

        model_outputs = dict(**inference_outputs, **generative_outputs)
        losses = self.loss(tensors, model_outputs, **loss_kwargs)

        return model_outputs, losses

    @abstractmethod
    def inference(
        self,
    ):
        pass

    @abstractmethod
    def generative(self):
        pass

    @abstractmethod
    def sample(self):
        pass
