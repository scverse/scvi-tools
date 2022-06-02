from typing import Any, Dict

import flax
import jax.numpy as jnp
from flax.core import FrozenDict
from flax.training import train_state
from jax import random

from scvi.utils._jax import device_selecting_PRNGKey

from ._base_module import JaxBaseModuleClass


class BatchTrainState(train_state.TrainState):
    batch_stats: FrozenDict[str, Any]


class JaxModuleWrapper:
    def __init__(
        self,
        module_cls: JaxBaseModuleClass,
        seed: int = 0,  # switch to using a global scvi.settings seed that gets forked everytime a modulewrapper is initialized by default
        use_gpu: bool = False,
        **module_kwargs,
    ) -> None:
        self.module_cls = module_cls
        self.module_kwargs = module_kwargs
        self._module = self._get_module()
        self._train_module = None
        self._eval_module = None
        self._train_state = None

        self.use_gpu = use_gpu
        self.key_fn = device_selecting_PRNGKey(use_gpu=self.use_gpu)
        self.seed_rng = self.key_fn(seed)
        self._set_rngs()

    @property
    def module(self):
        return self._module

    def _get_module(self, kwargs=None):
        kwargs = (
            self.module_kwargs if kwargs is None else {**self.module_kwargs, **kwargs}
        )
        return self.module_cls(**kwargs)

    def eval(self):
        if self._eval_module is None:
            self._eval_module = self._get_module(dict(is_training=False))
        self._module = self._eval_module

    def train(self):
        if self._train_module is None:
            self._train_module = self._get_module(dict(is_training=True))
        self._module = self._train_module

    @property
    def bound_module(self):
        return self.module.bind(
            {"params": self.params, "batch_stats": self.batch_stats},
            rngs=self.rngs,
        )

    # TODO(jhong): consider making this a mixin, or part of JaxBaseModule
    @property
    def apply(self):
        return self.module.apply

    @property
    def loss(self):
        return self.module.loss

    @property
    def init(self):
        return self.module.init

    @property
    def rngs(self) -> Dict[str, jnp.ndarray]:
        return self._split_rngs()

    def _set_rngs(self):
        required_rngs = self.module.required_rngs
        rng_keys = random.split(self.seed_rng, num=len(required_rngs) + 1)
        self.seed, module_rngs = rng_keys[0], rng_keys[1:]
        self._rngs = {k: module_rngs[i] for i, k in enumerate(required_rngs)}

    def _split_rngs(self):
        new_rngs = {}
        ret_rngs = {}
        for k, v in self._rngs.items():
            new_rngs[k], ret_rngs[k] = random.split(v)
        self._rngs = new_rngs
        return ret_rngs

    @property
    def train_state(self) -> BatchTrainState:
        return self._train_state

    @train_state.setter
    def train_state(self, train_state: BatchTrainState):
        self._train_state = train_state

    @property
    def params(self) -> flax.core.FrozenDict[str, Any]:
        return self.train_state.params

    @property
    def batch_stats(self) -> FrozenDict[str, Any]:
        return self.train_state.batch_stats


# should handle save and load
