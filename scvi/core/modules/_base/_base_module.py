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
        get_inference_input_kwargs={},
        get_generative_input_kwargs={},
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

        inference_inputs = self._get_inference_input(
            tensors, **get_inference_input_kwargs
        )
        inference_outputs = self.inference(**inference_inputs, **inference_kwargs)
        generative_inputs = self._get_generative_input(
            tensors, inference_outputs, **get_generative_input_kwargs
        )
        generative_outputs = self.generative(**generative_inputs, **generative_kwargs)

        model_outputs = dict(**inference_outputs, **generative_outputs)
        losses = self.loss(tensors, model_outputs, **loss_kwargs)

        return model_outputs, losses

    @abstractmethod
    def _get_inference_input(self, tensors):
        pass

    @abstractmethod
    def _get_generative_input(self, tensors, inference_outputs):
        pass

    @abstractmethod
    def inference(
        self,
        *args,
        **kwargs,
    ):
        pass

    @abstractmethod
    def generative(self, *args, **kwargs):
        pass

    @abstractmethod
    def loss(self, *args, **kwargs):
        pass

    @abstractmethod
    def sample(self):
        pass
