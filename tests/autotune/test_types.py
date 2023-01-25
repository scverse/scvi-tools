import scvi
from scvi.autotune import Tunable, TunableMixin


class DummyTrainingMixin:
    def train(self, lr: Tunable[int] = 1e-3):
        pass


class DummyDataSplitter(TunableMixin):
    def __init__(self, n_train: Tunable[int] = 1000, n_val: Tunable[int] = 100):
        self.n_train = n_train
        self.n_val = n_val


class DummyModel(TunableMixin, DummyTrainingMixin):

    _data_splitter_cls = DummyDataSplitter

    def __init__(
        self,
        n_input: int = 100,
        n_hidden: Tunable[int] = 128,
        n_latent: Tunable[int] = 10,
    ):
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.n_latent = n_latent


def test_tunable_mixin():
    model_cls = DummyModel
    manager = scvi.autotune.TunerManager(model_cls)

    # check all tunables are registered
    registry = manager._registry["tunables"]
    assert "n_train" in registry
    assert "n_val" in registry
    assert "n_hidden" in registry
    assert "n_latent" in registry
    assert "lr" in registry
    assert "n_input" not in registry
