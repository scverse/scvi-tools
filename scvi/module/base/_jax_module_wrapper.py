from typing import Any, Dict

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
        self.use_gpu = use_gpu
        self.key_fn = device_selecting_PRNGKey(use_gpu=self.use_gpu)
        self.seed_rng = self.key_fn(seed)
        self._set_rngs()

    def _get_module(self, kwargs=None):
        kwargs = (
            self.module_kwargs if kwargs is None else {**self.module_kwargs, **kwargs}
        )
        return self.module_cls(**kwargs)

    @property
    def module(self):
        if self._module is None:
            self._module = self._get_module()
        return self._module

    @property
    def rngs(self) -> Dict[str, jnp.ndarray]:
        return self._split_rngs()

    def _set_rngs(self):
        required_rngs = self.module_cls.required_rngs
        rng_keys = random.split(self.seed, num=len(required_rngs) + 1)
        self.seed, module_rngs = rng_keys[0], rng_keys[1:]
        self._rngs = {k: module_rngs[i] for i, k in enumerate(required_rngs)}

    def _split_rngs(self):
        new_rngs = {}
        ret_rngs = {}
        for k, v in self._rngs.items():
            new_rngs[k], ret_rngs[k] = random.split(v)
        self._rngs = new_rngs
        return ret_rngs


# should handle save and load
# should handle device management
# should handle train state
