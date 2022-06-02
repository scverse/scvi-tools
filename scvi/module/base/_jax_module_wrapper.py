from typing import Any

from flax.core import FrozenDict
from flax.training import train_state

from ._base_module import JaxBaseModuleClass
from scvi.utils._jax import device_selecting_PRNGKey


class BatchTrainState(train_state.TrainState):
    batch_stats: FrozenDict[str, Any]


class JaxModuleWrapper:
    def __init__(
        self,
        module_cls: JaxBaseModuleClass,
        seed: int = 0,  # switch to using a global scvi.settings seed that gets forked everytime a modulewrapper is initialized
        use_gpu: bool = False,
        **module_kwargs,
    ) -> None:
        self.module_cls = module_cls
        self.use_gpu = use_gpu
        self.seed = seed
        self.rng = device_selecting_PRNGKey(use_gpu=self.use_gpu)
        self.module_kwargs = module_kwargs

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


# should handle RNG params - get rid of set_rng
# should handle save and load
# should handle device management
# should handle train state
