from abc import abstractmethod
from typing import Tuple

import torch
import torch.nn as nn


class SCVILoss:
    def __init__(self, loss, reconstruction_loss, kl_local, kl_global):

        self._loss = loss if isinstance(loss, dict) else dict(loss=loss)
        self._reconstruction_loss = (
            reconstruction_loss
            if isinstance(reconstruction_loss, dict)
            else dict(reconstruction_loss=reconstruction_loss)
        )
        self._kl_local = (
            kl_local if isinstance(kl_local, dict) else dict(kl_local=kl_local)
        )
        self._kl_global = (
            kl_global if isinstance(kl_global, dict) else dict(kl_global=kl_global)
        )

    @staticmethod
    def _get_dict_sum(dictionary):
        sum = 0.0
        for value in dictionary.values():
            sum += value
        return sum

    @property
    def loss(self, sum=True):
        if sum:
            return self._get_dict_sum(self._loss)
        return self._loss

    @property
    def reconstruction_loss(self, sum=True):
        if sum:
            return self._get_dict_sum(self._reconstruction_loss)
        return self._reconstruction_loss

    @property
    def kl_local(self, sum=True):
        if sum:
            return self._get_dict_sum(self._kl_local)
        return self._kl_local

    @property
    def kl_global(self, sum=True):
        if sum:
            return self._get_dict_sum(self._kl_global)
        return self._kl_global

    @property
    def elbo(self):
        return


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
        losses = self.loss(
            tensors, inference_outputs, generative_outputs, **loss_kwargs
        )

        return inference_outputs, generative_outputs, losses

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
