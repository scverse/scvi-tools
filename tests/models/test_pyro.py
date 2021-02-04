import os

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
from pyro.infer.autoguide import AutoDiagonalNormal
from pyro.nn import PyroModule

from scvi import _CONSTANTS
from scvi.compose import PyroBaseModuleClass
from scvi.data import synthetic_iid
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import PyroTrainingPlan, Trainer


class BayesianRegressionPyroModel(PyroModule):
    def __init__(self, in_features, out_features):
        super().__init__()

        self.register_buffer("zero", torch.tensor(0.0))
        self.register_buffer("one", torch.tensor(1.0))
        self.register_buffer("ten", torch.tensor(10.0))

        self.linear = PyroModule[nn.Linear](in_features, out_features)

    def forward(self, x, y):
        sigma = pyro.sample("sigma", dist.Uniform(self.zero, self.ten))
        mean = self.linear(x).squeeze(-1)
        with pyro.plate("data", x.shape[0]):
            pyro.sample("obs", dist.Normal(mean, sigma), obs=y)
        return mean


class BayesianRegressionModule(PyroBaseModuleClass):
    def __init__(self, in_features, out_features):

        super().__init__()
        self.model = BayesianRegressionPyroModel(in_features, out_features)
        self.guide = AutoDiagonalNormal(self.model)

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x = tensor_dict[_CONSTANTS.X_KEY]
        y = tensor_dict[_CONSTANTS.LABELS_KEY]

        return (x, y.squeeze(1)), {}


def test_pyro_bayesian_regression(save_path):
    use_gpu = 0
    adata = synthetic_iid()
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(adata.shape[1], 1)
    plan = PyroTrainingPlan(model)
    trainer = Trainer(
        gpus=use_gpu,
        max_epochs=2,
    )
    trainer.fit(plan, train_dl)

    # test save and load
    post_dl = AnnDataLoader(adata, shuffle=False, batch_size=128)
    mean1 = []
    with torch.no_grad():
        for tensors in post_dl:
            args, kwargs = model._get_fn_args_from_batch(tensors)
            mean1.append(model(*args, **kwargs).cpu().numpy())
    mean1 = np.concatenate(mean1)

    model_save_path = os.path.join(save_path, "model_params.pt")
    torch.save(model.state_dict(), model_save_path)

    pyro.clear_param_store()
    new_model = BayesianRegressionModule(adata.shape[1], 1)
    # run model one step to get autoguide params
    try:
        new_model.load_state_dict(torch.load(model_save_path))
    except RuntimeError as err:
        if isinstance(new_model, PyroBaseModuleClass):
            plan = PyroTrainingPlan(new_model)
            trainer = Trainer(
                gpus=use_gpu,
                max_steps=1,
            )
            trainer.fit(plan, train_dl)
            new_model.load_state_dict(torch.load(model_save_path))
        else:
            raise err

    mean2 = []
    with torch.no_grad():
        for tensors in post_dl:
            args, kwargs = new_model._get_fn_args_from_batch(tensors)
            mean2.append(new_model(*args, **kwargs).cpu().numpy())
    mean2 = np.concatenate(mean2)

    np.testing.assert_array_equal(mean1, mean2)


def test_pyro_bayesian_regression_jit():
    use_gpu = 0
    adata = synthetic_iid()
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(adata.shape[1], 1)
    # warmup guide for JIT
    for tensors in train_dl:
        args, kwargs = model._get_fn_args_from_batch(tensors)
        model.guide(*args, **kwargs)
        break
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    plan = PyroTrainingPlan(model, loss_fn=pyro.infer.JitTrace_ELBO())
    trainer = Trainer(
        gpus=use_gpu,
        max_epochs=2,
    )
    trainer.fit(plan, train_dl)
