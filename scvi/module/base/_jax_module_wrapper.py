from typing import Any, Dict

import flax
import jax
import jax.numpy as jnp
from flax.core import FrozenDict
from flax.training import train_state
from jax import random
from jaxlib.xla_extension import Device

from scvi.utils._jax import device_selecting_PRNGKey

from ._base_module import JaxBaseModuleClass


class TrainStateWithState(train_state.TrainState):
    """TrainState with state attribute."""

    state: FrozenDict[str, Any]


class JaxModuleWrapper:
    """
    Wrapper class for Flax (Jax-backed) modules used to interact with model classes.

    This class maintains all state necessary for training and updating the state
    via the Flax module. The Flax module should remain stateless. In addition, the
    ``JaxModuleWrapper`` duck-types the methods of :class:`~scvi.module.base.BaseModuleClass`,
    which supports PyTorch-backed modules, to provide a consistent interface for
    :class:`~scvi.model.base.BaseModelClass`.

    Parameters
    ----------
    module_cls
        Flax module class to wrap.
    seed
        Random seed to initialize Jax RNGs with.
    **module_kwargs
        Keyword arguments that will be used to initialize ``module_cls``.
    """

    def __init__(
        self,
        module_cls: JaxBaseModuleClass,
        seed: int = 0,  # switch to using a global scvi.settings seed that gets forked everytime a modulewrapper is initialized by default
        **module_kwargs,
    ) -> None:
        self.module_cls = module_cls
        self.module_kwargs = module_kwargs
        self._module = self.module_cls(training=True, **self.module_kwargs)
        self._train_state = None

        self.key_fn = device_selecting_PRNGKey()
        self.seed_rng = self.key_fn(seed)
        self._set_rngs()

    @staticmethod
    def on_load(model):
        """
        Callback function run in :method:`~scvi.model.base.BaseModelClass.load` prior to loading module state dict.

        For some Pyro modules with AutoGuides, run one training step prior to loading state dict.
        """
        old_history = model.history_.copy()
        model.train(max_steps=1)
        model.history_ = old_history

    @property
    def device(self):  # noqa: D102
        return self.seed_rng.device()

    @property
    def module(self):  # noqa: D102
        return self._module

    def eval(self):
        """Switch to evaluation mode. Emulates Pytorch's interface."""
        self._module = self.module_cls(training=False, **self.module_kwargs)

    def train(self):
        """Switch to train mode. Emulates Pytorch's interface."""
        self._module = self.module_cls(training=True, **self.module_kwargs)

    @jax.jit
    def jit_inference(self, array_dict: Dict[str, jnp.ndarray], **kwargs):
        """
        Returns a method to run base module's inference method

        Parameters
        ----------
        **kwargs
            Inference kwargs
        """
        vars_in = {"params": self.params, **self.state}
        inference_input = self.module._get_inference_input(array_dict)
        out = self.module.apply(
            vars_in, rngs=self.rngs, method=self.module.inference, **inference_input
        )
        return out

    @jax.jit
    def jit_generative(self, array_dict: Dict[str, jnp.ndarray], **kwargs):
        """
        Returns a method to run base module's generative method

        Parameters
        ----------
        **kwargs
            generative kwargs
        """
        vars_in = {"params": self.params, **self.state}
        generative_input = self.module._get_generative_input(array_dict)
        out = self.module.apply(
            vars_in, rngs=self.rngs, method=self.module.generative, **generative_input
        )
        return out

    @jax.jit
    def jit_forward(self, array_dict: Dict[str, jnp.ndarray], **kwargs):
        """
        Returns a method to run base module's forward (__call__) method

        Parameters
        ----------
        **kwargs
            Inference kwargs
        """
        vars_in = {"params": self.params, **self.state}
        out = self.module.apply(vars_in, array_dict, rngs=self.rngs, **kwargs)
        return out

    @property
    def apply(self):
        """Apply function of the Flax module."""
        return self.module.apply

    @property
    def loss(self):
        """Loss function of the Flax module."""
        return self.module.loss

    @property
    def init(self):
        """Init function of the Flax module."""
        return self.module.init

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
        required_rngs = self.module.required_rngs
        rng_keys = random.split(self.seed_rng, num=len(required_rngs) + 1)
        self.seed, module_rngs = rng_keys[0], rng_keys[1:]
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
    def train_state(self) -> TrainStateWithState:
        """Train state containing learned parameter values from training."""
        return self._train_state

    @train_state.setter
    def train_state(self, train_state: TrainStateWithState):
        self._train_state = train_state

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
        self._check_train_state_is_not_none()
        self.train_state = flax.serialization.from_state_dict(
            self.train_state, state_dict
        )

    def to(self, device: Device):
        """Move module to device."""
        self._check_train_state_is_not_none()
        self.train_state = jax.tree_util.tree_map(
            lambda x: jax.device_put(x, device), self.train_state
        )
        self.seed_rng = jax.device_put(self.seed_rng, device)
        self._rngs = jax.device_put(self._rngs, device)

    def _check_train_state_is_not_none(self):
        if self.train_state is None:
            raise RuntimeError(
                "Train state is not set. Module has not been trained. If loading, initialize the model."
            )
