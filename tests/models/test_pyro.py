import os

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
from pyro.infer.autoguide import AutoDiagonalNormal
from pyro.nn import PyroModule, PyroSample
from pytorch_lightning.callbacks import Callback

from scvi import _CONSTANTS
from scvi.compose import PyroBaseModuleClass
from scvi.data import synthetic_iid
from scvi.dataloaders import AnnDataLoader
from scvi.lightning import PyroTrainingPlan, Trainer


class PyroSampleCallback(Callback):
    def on_train_start(self, trainer, pl_module):
        """
        Way to use PyroSample and have it be on GPU.

        In this case, the PyroSample objects are added after the
        model has already been moved to cuda.
        """

        pyro_model = pl_module.module.model
        pyro_model.linear.weight = PyroSample(
            dist.Normal(pyro_model.zero, pyro_model.one)
            .expand([pyro_model.out_features, pyro_model.in_features])
            .to_event(2)
        )
        pyro_model.linear.bias = PyroSample(
            dist.Normal(pyro_model.zero, pyro_model.ten)
            .expand([pyro_model.out_features])
            .to_event(1)
        )


class BayesianRegressionPyroModel(PyroModule):
    def __init__(self, in_features, out_features):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features

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
    use_gpu = int(torch.cuda.is_available())
    adata = synthetic_iid()
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(adata.shape[1], 1)
    plan = PyroTrainingPlan(model)
    trainer = Trainer(gpus=use_gpu, max_epochs=2, callbacks=[PyroSampleCallback()])
    trainer.fit(plan, train_dl)

    # test save and load
    quants = model.guide.quantiles([0.5])
    sigma_median = quants["sigma"][0].cpu().detach().numpy()
    linear_median = quants["linear.weight"][0].cpu().detach().numpy()

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
                gpus=use_gpu, max_steps=1, callbacks=[PyroSampleCallback()]
            )
            trainer.fit(plan, train_dl)
            new_model.load_state_dict(torch.load(model_save_path))
        else:
            raise err

    quants = new_model.guide.quantiles([0.5])
    sigma_median_new = quants["sigma"][0].cpu().detach().numpy()
    linear_median_new = quants["linear.weight"][0].cpu().detach().numpy()

    np.testing.assert_array_equal(sigma_median_new, sigma_median)
    np.testing.assert_array_equal(linear_median_new, linear_median)


def test_pyro_bayesian_regression_jit():
    use_gpu = int(torch.cuda.is_available())
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
