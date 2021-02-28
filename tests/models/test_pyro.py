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
from scvi.train import PyroTrainingPlan, Trainer


class PyroJitGuideWarmup(Callback):
    def __init__(self, train_dl) -> None:
        super().__init__()
        self.dl = train_dl

    def on_train_start(self, trainer, pl_module):
        """
        Way to warmup Pyro Guide in an automated way.

        Also device agnostic.
        """

        # warmup guide for JIT
        pyro_model = pl_module.module.model
        dev = pyro_model.linear.weight.device
        pyro_guide = pl_module.module.guide
        for tensors in self.dl:
            tens = {k: t.to(dev) for k, t in tensors.items()}
            args, kwargs = pl_module.module._get_fn_args_from_batch(tens)
            pyro_guide(*args, **kwargs)
            break


class BayesianRegressionPyroModel(PyroModule):
    def __init__(self, in_features, out_features):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features

        self.register_buffer("zero", torch.tensor(0.0))
        self.register_buffer("one", torch.tensor(1.0))
        self.register_buffer("ten", torch.tensor(10.0))

        self.linear = PyroModule[nn.Linear](in_features, out_features)
        self.linear.weight = PyroSample(
            lambda prior: dist.Normal(self.zero, self.one)
            .expand([self.out_features, self.in_features])
            .to_event(2)
        )
        self.linear.bias = PyroSample(
            lambda prior: dist.Normal(self.zero, self.ten)
            .expand([self.out_features])
            .to_event(1)
        )

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
    trainer = Trainer(
        gpus=use_gpu,
        max_epochs=2,
    )
    trainer.fit(plan, train_dl)
    if use_gpu == 1:
        model.cuda()

    # test Predictive
    num_samples = 5
    predictive = model.create_predictive(num_samples=num_samples)
    for tensor_dict in train_dl:
        args, kwargs = model._get_fn_args_from_batch(tensor_dict)
        _ = {
            k: v.detach().cpu().numpy()
            for k, v in predictive(*args, **kwargs).items()
            if k != "obs"
        }
    # test save and load
    # cpu/gpu has minor difference
    model.cpu()
    quants = model.guide.quantiles([0.5])
    sigma_median = quants["sigma"][0].detach().cpu().numpy()
    linear_median = quants["linear.weight"][0].detach().cpu().numpy()

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

    quants = new_model.guide.quantiles([0.5])
    sigma_median_new = quants["sigma"][0].detach().cpu().numpy()
    linear_median_new = quants["linear.weight"][0].detach().cpu().numpy()

    np.testing.assert_array_equal(sigma_median_new, sigma_median)
    np.testing.assert_array_equal(linear_median_new, linear_median)


def test_pyro_bayesian_regression_jit():
    use_gpu = int(torch.cuda.is_available())
    adata = synthetic_iid()
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(adata.shape[1], 1)
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    plan = PyroTrainingPlan(model, loss_fn=pyro.infer.JitTrace_ELBO())
    trainer = Trainer(
        gpus=use_gpu, max_epochs=2, callbacks=[PyroJitGuideWarmup(train_dl)]
    )
    trainer.fit(plan, train_dl)

    # 100 features, 1 for sigma, 1 for bias
    assert list(model.guide.parameters())[0].shape[0] == 102

    if use_gpu == 1:
        model.cuda()

    # test Predictive
    num_samples = 5
    predictive = model.create_predictive(num_samples=num_samples)
    for tensor_dict in train_dl:
        args, kwargs = model._get_fn_args_from_batch(tensor_dict)
        _ = {
            k: v.detach().cpu().numpy()
            for k, v in predictive(*args, **kwargs).items()
            if k != "obs"
        }
