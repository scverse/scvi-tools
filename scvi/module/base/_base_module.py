from __future__ import annotations

import warnings
from abc import abstractmethod
from dataclasses import field
from typing import Any, Callable, Dict, Iterable, Optional, Tuple, Union

import chex
import flax
import jax
import jax.numpy as jnp
import numpy as np
import pyro
import torch
from flax.core import FrozenDict
from flax.training import train_state
from jax import random
from jaxlib.xla_extension import Device
from numpyro.distributions import Distribution
from pyro.infer.predictive import Predictive
from torch import nn

from scvi import settings
from scvi._types import LatentDataType, LossRecord, Tensor
from scvi.autotune._types import TunableMixin
from scvi.utils._jax import device_selecting_PRNGKey

from ._decorators import auto_move_data
from ._pyro import AutoMoveDataPredictive


class LossRecorder:
    """
    Loss signature for models.

    This class provides an organized way to record the model loss, as well as
    the components of the ELBO. This may also be used in MLE, MAP, EM methods.
    The loss is used for backpropagation during inference. The other parameters
    are used for logging/early stopping during inference.

    Parameters
    ----------
    loss
        Tensor with loss for minibatch. Should be one dimensional with one value.
        Note that loss should be a :class:`~torch.Tensor` and not the result of ``.item()``.
    reconstruction_loss
        Reconstruction loss for each observation in the minibatch. If a tensor, converted to
        a dictionary with key "reconstruction_loss" and value as tensor
    kl_local
        KL divergence associated with each observation in the minibatch. If a tensor, converted to
        a dictionary with key "kl_local" and value as tensor
    kl_global
        Global kl divergence term. Should be one dimensional with one value. If a tensor, converted to
        a dictionary with key "kl_global" and value as tensor
    **kwargs
        Additional metrics can be passed as keyword arguments and will
        be available as attributes of the object.
    """

    def __init__(
        self,
        loss: LossRecord,
        reconstruction_loss: Optional[LossRecord] = None,
        kl_local: Optional[LossRecord] = None,
        kl_global: Optional[LossRecord] = None,
        **kwargs,
    ):
        warnings.warn(
            "LossRecorder is deprecated and will be removed in version 0.20.0. Please use LossOutput",
            category=DeprecationWarning,
        )
        self._loss_output = LossOutput(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            kl_local=kl_local,
            kl_global=kl_global,
            extra_metrics=kwargs,
        )
        self.extra_metric_attrs = []
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.extra_metric_attrs.append(key)

    @property
    def loss(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self._loss_output.loss

    @property
    def reconstruction_loss(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self.dict_sum(self._loss_output.reconstruction_loss)

    @property
    def _reconstruction_loss(self):
        return self._loss_output.reconstruction_loss

    @property
    def kl_local(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self.dict_sum(self._loss_output.kl_local)

    @property
    def _kl_local(self):
        return self._loss_output.kl_local

    @property
    def reconstruction_loss_sum(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self._loss_output.reconstruction_loss_sum

    @property
    def kl_local_sum(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self._loss_output.kl_local_sum

    @property
    def kl_global_sum(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self._loss_output.kl_global_sum

    @property
    def kl_global(self) -> Union[torch.Tensor, jnp.ndarray]:  # noqa: D102
        return self.dict_sum(self._loss_output.kl_global)

    def dict_sum(self, x):
        """Wrapper of LossOutput.dict_sum."""
        return self._loss_output.dict_sum(x)


@chex.dataclass
class LossOutput:
    """
    Loss signature for models.

    This class provides an organized way to record the model loss, as well as
    the components of the ELBO. This may also be used in MLE, MAP, EM methods.
    The loss is used for backpropagation during inference. The other parameters
    are used for logging/early stopping during inference.

    Parameters
    ----------
    loss
        Tensor with loss for minibatch. Should be one dimensional with one value.
        Note that loss should be in an array/tensor and not a float.
    reconstruction_loss
        Reconstruction loss for each observation in the minibatch. If a tensor, converted to
        a dictionary with key "reconstruction_loss" and value as tensor.
    kl_local
        KL divergence associated with each observation in the minibatch. If a tensor, converted to
        a dictionary with key "kl_local" and value as tensor.
    kl_global
        Global KL divergence term. Should be one dimensional with one value. If a tensor, converted to
        a dictionary with key "kl_global" and value as tensor.
    extra_metrics
        Additional metrics can be passed as arrays/tensors or dictionaries of
        arrays/tensors.
    n_obs_minibatch
        Number of observations in the minibatch. If None, will be inferred from
        the shape of the reconstruction_loss tensor.


    Examples
    --------
    >>> loss_output = LossOutput(
    ...     loss=loss,
    ...     reconstruction_loss=reconstruction_loss,
    ...     kl_local=kl_local,
    ...     extra_metrics={"x": scalar_tensor_x, "y": scalar_tensor_y},
    ... )
    """

    loss: LossRecord
    reconstruction_loss: Optional[LossRecord] = None
    kl_local: Optional[LossRecord] = None
    kl_global: Optional[LossRecord] = None
    extra_metrics: Optional[Dict[str, Tensor]] = field(default_factory=dict)
    n_obs_minibatch: Optional[int] = None
    reconstruction_loss_sum: Tensor = field(default=None, init=False)
    kl_local_sum: Tensor = field(default=None, init=False)
    kl_global_sum: Tensor = field(default=None, init=False)

    def __post_init__(self):
        self.loss = self.dict_sum(self.loss)

        if self.n_obs_minibatch is None and self.reconstruction_loss is None:
            raise ValueError(
                "Must provide either n_obs_minibatch or reconstruction_loss"
            )

        default = 0 * self.loss
        if self.reconstruction_loss is None:
            self.reconstruction_loss = default
        if self.kl_local is None:
            self.kl_local = default
        if self.kl_global is None:
            self.kl_global = default
        self.reconstruction_loss = self._as_dict("reconstruction_loss")
        self.kl_local = self._as_dict("kl_local")
        self.kl_global = self._as_dict("kl_global")
        self.reconstruction_loss_sum = self.dict_sum(self.reconstruction_loss).sum()
        self.kl_local_sum = self.dict_sum(self.kl_local).sum()
        self.kl_global_sum = self.dict_sum(self.kl_global)

        if self.reconstruction_loss is not None and self.n_obs_minibatch is None:
            rec_loss = self.reconstruction_loss
            self.n_obs_minibatch = list(rec_loss.values())[0].shape[0]

    @staticmethod
    def dict_sum(dictionary: Union[Dict[str, Tensor], Tensor]):
        """Sum over elements of a dictionary."""
        if isinstance(dictionary, dict):
            return sum(dictionary.values())
        else:
            return dictionary

    @property
    def extra_metrics_keys(self) -> Iterable[str]:
        """Keys for extra metrics."""
        return self.extra_metrics.keys()

    def _as_dict(self, attr_name: str):
        attr = getattr(self, attr_name)
        if isinstance(attr, dict):
            return attr
        else:
            return {attr_name: attr}


class BaseModuleClass(TunableMixin, nn.Module):
    """Abstract class for scvi-tools modules."""

    def __init__(
        self,
    ):
        super().__init__()

    @property
    def device(self):  # noqa: D102
        device = list({p.device for p in self.parameters()})
        if len(device) > 1:
            raise RuntimeError("Module tensors on multiple devices.")
        return device[0]

    def on_load(self, model):
        """Callback function run in :meth:`~scvi.model.base.BaseModelClass.load` prior to loading module state dict."""

    @auto_move_data
    def forward(
        self,
        tensors,
        get_inference_input_kwargs: Optional[dict] = None,
        get_generative_input_kwargs: Optional[dict] = None,
        inference_kwargs: Optional[dict] = None,
        generative_kwargs: Optional[dict] = None,
        loss_kwargs: Optional[dict] = None,
        compute_loss=True,
    ) -> Union[
        Tuple[torch.Tensor, torch.Tensor],
        Tuple[torch.Tensor, torch.Tensor, LossRecorder],
    ]:
        """
        Forward pass through the network.

        Parameters
        ----------
        tensors
            tensors to pass through
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
            Whether to compute loss on forward pass. This adds
            another return value.
        """
        return _generic_forward(
            self,
            tensors,
            inference_kwargs,
            generative_kwargs,
            loss_kwargs,
            get_inference_input_kwargs,
            get_generative_input_kwargs,
            compute_loss,
        )

    @abstractmethod
    def _get_inference_input(self, tensors: Dict[str, torch.Tensor], **kwargs):
        """Parse tensors dictionary for inference related values."""

    @abstractmethod
    def _get_generative_input(
        self,
        tensors: Dict[str, torch.Tensor],
        inference_outputs: Dict[str, torch.Tensor],
        **kwargs,
    ):
        """Parse tensors dictionary for generative related values."""

    @abstractmethod
    def inference(
        self,
        *args,
        **kwargs,
    ) -> Dict[str, Union[torch.Tensor, torch.distributions.Distribution]]:
        """
        Run the recognition model.

        In the case of variational inference, this function will perform steps related to
        computing variational distribution parameters. In a VAE, this will involve running
        data through encoder networks.

        This function should return a dictionary with str keys and :class:`~torch.Tensor` values.
        """

    @abstractmethod
    def generative(
        self, *args, **kwargs
    ) -> Dict[str, Union[torch.Tensor, torch.distributions.Distribution]]:
        """
        Run the generative model.

        This function should return the parameters associated with the likelihood of the data.
        This is typically written as :math:`p(x|z)`.

        This function should return a dictionary with str keys and :class:`~torch.Tensor` values.
        """

    @abstractmethod
    def loss(self, *args, **kwargs) -> LossOutput:
        """
        Compute the loss for a minibatch of data.

        This function uses the outputs of the inference and generative functions to compute
        a loss. This many optionally include other penalty terms, which should be computed here.

        This function should return an object of type :class:`~scvi.module.base.LossOutput`.
        """

    @abstractmethod
    def sample(self, *args, **kwargs):
        """Generate samples from the learned model."""


class BaseLatentModeModuleClass(BaseModuleClass):
    """Abstract base class for scvi-tools modules that support latent mode."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._latent_data_type = None

    @property
    def latent_data_type(self) -> Union[LatentDataType, None]:
        """The latent data type associated with this module."""
        return self._latent_data_type

    @latent_data_type.setter
    def latent_data_type(self, latent_data_type):
        """Set latent data type associated with this module."""
        self._latent_data_type = latent_data_type

    @abstractmethod
    def _cached_inference(self, *args, **kwargs):
        """Uses the cached latent mode distribution to perform inference, thus bypassing the encoder."""

    @abstractmethod
    def _regular_inference(self, *args, **kwargs):
        """Runs inference (encoder forward pass)."""

    @auto_move_data
    def inference(self, *args, **kwargs):
        """
        Main inference call site.

        Branches off to regular or cached inference depending on the latent data
        type of the module.
        """
        if self.latent_data_type is None:
            return self._regular_inference(*args, **kwargs)
        else:
            return self._cached_inference(*args, **kwargs)


def _get_dict_if_none(param):
    param = {} if not isinstance(param, dict) else param

    return param


class PyroBaseModuleClass(TunableMixin, nn.Module):
    """
    Base module class for Pyro models.

    In Pyro, ``model`` and ``guide`` should have the same signature. Out of convenience,
    the forward function of this class passes through to the forward of the ``model``.

    There are two ways this class can be equipped with a model and a guide. First,
    ``model`` and ``guide`` can be class attributes that are :class:`~pyro.nn.PyroModule`
    instances. The implemented ``model`` and ``guide`` class method can then return the (private) attributes.
    Second, ``model`` and ``guide`` methods can be written directly (see Pyro scANVI example)
    https://pyro.ai/examples/scanvi.html.

    The ``model`` and ``guide`` may also be equipped with ``n_obs`` attributes, which can be set
    to ``None`` (e.g., ``self.n_obs = None``). This attribute may be helpful in designating the
    size of observation-specific Pyro plates. The value will be updated automatically by
    :class:`~scvi.train.PyroTrainingPlan`, provided that it is given the number of training examples
    upon initialization.

    Parameters
    ----------
    on_load_kwargs
        Dictionary containing keyword args to use in ``self.on_load``.
    """

    def __init__(self, on_load_kwargs: Optional[dict] = None):
        super().__init__()
        self.on_load_kwargs = on_load_kwargs or {}

    @staticmethod
    @abstractmethod
    def _get_fn_args_from_batch(
        tensor_dict: Dict[str, torch.Tensor]
    ) -> Union[Iterable, dict]:
        """
        Parse the minibatched data to get the correct inputs for ``model`` and ``guide``.

        In Pyro, ``model`` and ``guide`` must have the same signature. This is a helper method
        that gets the args and kwargs for these two methods. This helper method aids ``forward`` and
        ``guide`` in having transparent signatures, as well as allows use of our generic
        :class:`~scvi.dataloaders.AnnDataLoader`.

        Returns
        -------
        args and kwargs for the functions, args should be an Iterable and kwargs a dictionary.
        """

    @property
    @abstractmethod
    def model(self):  # noqa: D102
        pass

    @property
    @abstractmethod
    def guide(self):  # noqa: D102
        pass

    @property
    def list_obs_plate_vars(self):
        """
        Model annotation for minibatch training with pyro plate.

        A dictionary with:
        1. "name" - the name of observation/minibatch plate;
        2. "in" - indexes of model args to provide to encoder network when using amortised inference;
        3. "sites" - dictionary with
            keys - names of variables that belong to the observation plate (used to recognise
             and merge posterior samples for minibatch variables)
            values - the dimensions in non-plate axis of each variable (used to construct output
             layer of encoder network when using amortised inference)
        """
        return {"name": "", "in": [], "sites": {}}

    def on_load(self, model):
        """
        Callback function run in :method:`~scvi.model.base.BaseModelClass.load` prior to loading module state dict.

        For some Pyro modules with AutoGuides, run one training step prior to loading state dict.
        """
        old_history = model.history_.copy()
        model.train(max_steps=1, **self.on_load_kwargs)
        model.history_ = old_history
        pyro.clear_param_store()

    def create_predictive(
        self,
        model: Optional[Callable] = None,
        posterior_samples: Optional[dict] = None,
        guide: Optional[Callable] = None,
        num_samples: Optional[int] = None,
        return_sites: Tuple[str] = (),
        parallel: bool = False,
    ) -> Predictive:
        """
        Creates a :class:`~pyro.infer.Predictive` object.

        Parameters
        ----------
        model
            Python callable containing Pyro primitives. Defaults to ``self.model``.
        posterior_samples
            Dictionary of samples from the posterior
        guide
            Optional guide to get posterior samples of sites not present
            in ``posterior_samples``. Defaults to ``self.guide``
        num_samples
            Number of samples to draw from the predictive distribution.
            This argument has no effect if ``posterior_samples`` is non-empty, in which case,
            the leading dimension size of samples in ``posterior_samples`` is used.
        return_sites
            Sites to return; by default only sample sites not present
            in ``posterior_samples`` are returned.
        parallel
            predict in parallel by wrapping the existing model
            in an outermost ``plate`` messenger. Note that this requires that the model has
            all batch dims correctly annotated via :class:`~pyro.plate`.
        """
        if model is None:
            model = self.model
        if guide is None:
            guide = self.guide
        predictive = AutoMoveDataPredictive(
            model=model,
            posterior_samples=posterior_samples,
            guide=guide,
            num_samples=num_samples,
            return_sites=return_sites,
            parallel=parallel,
        )
        # necessary to comply with auto_move_data decorator
        predictive.eval()

        return predictive

    def forward(self, *args, **kwargs):
        """Passthrough to Pyro model."""
        return self.model(*args, **kwargs)


class TrainStateWithState(train_state.TrainState):
    """TrainState with state attribute."""

    state: FrozenDict[str, Any]


class JaxBaseModuleClass(TunableMixin, flax.linen.Module):
    """
    Abstract class for Jax-based scvi-tools modules.

    The :class:`~scvi.module.base.JaxBaseModuleClass` provides an interface for Jax-backed
    modules consistent with the :class:`~scvi.module.base.BaseModuleClass`.

    Any subclass must has a `training` parameter in its constructor, as well as
    use the `@flax_configure` decorator.

    Children of :class:`~scvi.module.base.JaxBaseModuleClass` should
    use the instance attribute ``self.training`` to appropriately modify
    the behavior of the model whether it is in training or evaluation mode.
    """

    def configure(self) -> None:
        """Add necessary attrs."""
        self.training = None
        self.train_state = None
        self.seed = settings.seed
        self.seed_rng = device_selecting_PRNGKey()(self.seed)
        self._set_rngs()

    @abstractmethod
    def setup(self):
        """
        Flax setup method.

        With scvi-tools we prefer to use the setup parameterization of
        flax.linen Modules. This lends the interface to be more like
        PyTorch. More about this can be found here:

        https://flax.readthedocs.io/en/latest/design_notes/setup_or_nncompact.html
        """

    @property
    @abstractmethod
    def required_rngs(self):
        """Returns a tuple of rng sequence names required for this Flax module."""
        return ("params",)

    def __call__(
        self,
        tensors: Dict[str, jnp.ndarray],
        get_inference_input_kwargs: Optional[dict] = None,
        get_generative_input_kwargs: Optional[dict] = None,
        inference_kwargs: Optional[dict] = None,
        generative_kwargs: Optional[dict] = None,
        loss_kwargs: Optional[dict] = None,
        compute_loss=True,
    ) -> Union[
        Tuple[jnp.ndarray, jnp.ndarray],
        Tuple[jnp.ndarray, jnp.ndarray, LossRecorder],
    ]:
        """
        Forward pass through the network.

        Parameters
        ----------
        tensors
            tensors to pass through
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
            Whether to compute loss on forward pass. This adds
            another return value.
        """
        return _generic_forward(
            self,
            tensors,
            inference_kwargs,
            generative_kwargs,
            loss_kwargs,
            get_inference_input_kwargs,
            get_generative_input_kwargs,
            compute_loss,
        )

    @abstractmethod
    def _get_inference_input(self, tensors: Dict[str, jnp.ndarray], **kwargs):
        """Parse tensors dictionary for inference related values."""

    @abstractmethod
    def _get_generative_input(
        self,
        tensors: Dict[str, jnp.ndarray],
        inference_outputs: Dict[str, jnp.ndarray],
        **kwargs,
    ):
        """Parse tensors dictionary for generative related values."""

    @abstractmethod
    def inference(
        self,
        *args,
        **kwargs,
    ) -> Dict[str, Union[jnp.ndarray, Distribution]]:
        """
        Run the recognition model.

        In the case of variational inference, this function will perform steps related to
        computing variational distribution parameters. In a VAE, this will involve running
        data through encoder networks.

        This function should return a dictionary with str keys and :class:`~jnp.ndarray` values.
        """

    @abstractmethod
    def generative(
        self, *args, **kwargs
    ) -> Dict[str, Union[jnp.ndarray, Distribution]]:
        """
        Run the generative model.

        This function should return the parameters associated with the likelihood of the data.
        This is typically written as :math:`p(x|z)`.

        This function should return a dictionary with str keys and :class:`~jnp.ndarray` values.
        """

    @abstractmethod
    def loss(self, *args, **kwargs) -> LossOutput:
        """
        Compute the loss for a minibatch of data.

        This function uses the outputs of the inference and generative functions to compute
        a loss. This many optionally include other penalty terms, which should be computed here.

        This function should return an object of type :class:`~scvi.module.base.LossOutput`.
        """

    @property
    def device(self):  # noqa: D102
        return self.seed_rng.device()

    def train(self):
        """Switch to train mode. Emulates Pytorch's interface."""
        self.training = True

    def eval(self):
        """Switch to evaluation mode. Emulates Pytorch's interface."""
        self.training = False

    @property
    def rngs(self) -> Dict[str, jnp.ndarray]:
        """
        Dictionary of RNGs mapping required RNG name to RNG values.

        Calls ``self._split_rngs()`` resulting in newly generated RNGs on
        every reference to ``self.rngs``.
        """
        return self._split_rngs()

    def _set_rngs(self):
        """Creates RNGs split off of the seed RNG for each RNG required by the module."""
        required_rngs = self.required_rngs
        rng_keys = random.split(self.seed_rng, num=len(required_rngs) + 1)
        self.seed_rng, module_rngs = rng_keys[0], rng_keys[1:]
        self._rngs = {k: module_rngs[i] for i, k in enumerate(required_rngs)}

    def _split_rngs(self):
        """
        Regenerates the current set of RNGs and returns newly split RNGs.

        Importantly, this method does not reuse RNGs in future references to ``self.rngs``.
        """
        new_rngs = {}
        ret_rngs = {}
        for k, v in self._rngs.items():
            new_rngs[k], ret_rngs[k] = random.split(v)
        self._rngs = new_rngs
        return ret_rngs

    @property
    def params(self) -> FrozenDict[str, Any]:  # noqa: D102
        self._check_train_state_is_not_none()
        return self.train_state.params

    @property
    def state(self) -> FrozenDict[str, Any]:  # noqa: D102
        self._check_train_state_is_not_none()
        return self.train_state.state

    def state_dict(self) -> Dict[str, Any]:
        """Returns a serialized version of the train state as a dictionary."""
        self._check_train_state_is_not_none()
        return flax.serialization.to_state_dict(self.train_state)

    def load_state_dict(self, state_dict: Dict[str, Any]):
        """Load a state dictionary into a train state."""
        if self.train_state is None:
            raise RuntimeError(
                "Train state is not set. Train for one iteration prior to loading state dict."
            )
        self.train_state = flax.serialization.from_state_dict(
            self.train_state, state_dict
        )

    def to(self, device: Device):
        """Move module to device."""
        if device is not self.device:
            if self.train_state is not None:
                self.train_state = jax.tree_util.tree_map(
                    lambda x: jax.device_put(x, device), self.train_state
                )

            self.seed_rng = jax.device_put(self.seed_rng, device)
            self._rngs = jax.device_put(self._rngs, device)

    def _check_train_state_is_not_none(self):
        if self.train_state is None:
            raise RuntimeError("Train state is not set. Module has not been trained.")

    def as_bound(self) -> JaxBaseModuleClass:
        """Module bound with parameters learned from training."""
        return self.bind(
            {"params": self.params, **self.state},
            rngs=self.rngs,
        )

    def get_jit_inference_fn(
        self,
        get_inference_input_kwargs: Optional[Dict[str, Any]] = None,
        inference_kwargs: Optional[Dict[str, Any]] = None,
    ) -> Callable[
        [Dict[str, jnp.ndarray], Dict[str, jnp.ndarray]], Dict[str, jnp.ndarray]
    ]:
        """
        Create a method to run inference using the bound module.

        Parameters
        ----------
        get_inference_input_kwargs
            Keyword arguments to pass to subclass `_get_inference_input`
        inference_kwargs
            Keyword arguments  for subclass `inference` method

        Returns
        -------
        A callable taking rngs and array_dict as input and returning the output
        of the `inference` method. This callable runs `_get_inference_input`.
        """
        vars_in = {"params": self.params, **self.state}
        get_inference_input_kwargs = _get_dict_if_none(get_inference_input_kwargs)
        inference_kwargs = _get_dict_if_none(inference_kwargs)

        @jax.jit
        def _run_inference(rngs, array_dict):
            module = self.clone()
            inference_input = module._get_inference_input(array_dict)
            out = module.apply(
                vars_in,
                rngs=rngs,
                method=module.inference,
                **inference_input,
                **inference_kwargs,
            )
            return out

        return _run_inference

    @staticmethod
    def on_load(model):
        """
        Callback function run in :meth:`~scvi.model.base.BaseModelClass.load` prior to loading module state dict.

        Run one training step prior to loading state dict in order to initialize params.
        """
        old_history = model.history_.copy()
        model.train(max_steps=1)
        model.history_ = old_history

    @staticmethod
    def as_numpy_array(x: jnp.ndarray):
        """Converts a jax device array to a numpy array."""
        return np.array(jax.device_get(x))


def _generic_forward(
    module,
    tensors,
    inference_kwargs,
    generative_kwargs,
    loss_kwargs,
    get_inference_input_kwargs,
    get_generative_input_kwargs,
    compute_loss,
):
    """Core of the forward call shared by PyTorch- and Jax-based modules."""
    inference_kwargs = _get_dict_if_none(inference_kwargs)
    generative_kwargs = _get_dict_if_none(generative_kwargs)
    loss_kwargs = _get_dict_if_none(loss_kwargs)
    get_inference_input_kwargs = _get_dict_if_none(get_inference_input_kwargs)
    get_generative_input_kwargs = _get_dict_if_none(get_generative_input_kwargs)

    inference_inputs = module._get_inference_input(
        tensors, **get_inference_input_kwargs
    )
    inference_outputs = module.inference(**inference_inputs, **inference_kwargs)
    generative_inputs = module._get_generative_input(
        tensors, inference_outputs, **get_generative_input_kwargs
    )
    generative_outputs = module.generative(**generative_inputs, **generative_kwargs)
    if compute_loss:
        losses = module.loss(
            tensors, inference_outputs, generative_outputs, **loss_kwargs
        )
        return inference_outputs, generative_outputs, losses
    else:
        return inference_outputs, generative_outputs
