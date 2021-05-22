import os

import matplotlib.pyplot as plt
import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn as nn
from anndata import AnnData
from pyro import clear_param_store
from pyro.infer.autoguide import AutoNormal, init_to_mean
from pyro.nn import PyroModule, PyroSample

from scvi import _CONSTANTS
from scvi.data import register_tensor_from_anndata, synthetic_iid
from scvi.dataloaders import AnnDataLoader
from scvi.model.base import (
    BaseModelClass,
    PyroJitGuideWarmup,
    PyroSampleMixin,
    PyroSviTrainMixin,
)
from scvi.module.base import PyroBaseModuleClass
from scvi.train import PyroTrainingPlan, Trainer


class BayesianRegressionPyroModel(PyroModule):
    def __init__(self, in_features, out_features, per_cell_weight=False):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.n_obs = None

        self.per_cell_weight = per_cell_weight

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

    def create_plates(self, x, y, ind_x):
        """Function for creating plates is needed when using AutoGuides
        and should have the same call signature as model"""
        return pyro.plate("obs_plate", size=self.n_obs, dim=-2, subsample=ind_x)

    def list_obs_plate_vars(self):
        """Model annotation for minibatch training with pyro plate.

        A dictionary with:
        1. "name" - the name of observation/minibatch plate;
        2. "in" - indexes of model args to provide to encoder network when using amortised inference;
        3. "sites" - dictionary with
            keys - names of variables that belong to the observation plate (used to recognise
             and merge posterior samples for minibatch variables)
            values - the dimensions in non-plate axis of each variable (used to construct output
             layer of encoder network when using amortised inference)
        """
        return {
            "name": "obs_plate",
            "in": [0],  # model args index for expression data
            "sites": {"per_cell_weights": 1},
        }

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x = tensor_dict[_CONSTANTS.X_KEY]
        y = tensor_dict[_CONSTANTS.LABELS_KEY]
        ind_x = tensor_dict["ind_x"].long().squeeze()
        return (x, y, ind_x), {}

    def forward(self, x, y, ind_x):

        obs_plate = self.create_plates(x, y, ind_x)

        sigma = pyro.sample("sigma", dist.Exponential(self.one))

        mean = self.linear(x).squeeze(-1)

        if self.per_cell_weight:
            with obs_plate:
                per_cell_weights = pyro.sample(
                    "per_cell_weights", dist.Normal(self.zero, self.one)
                )
                mean = mean + per_cell_weights.squeeze(-1)

        with obs_plate:
            pyro.sample("obs", dist.Normal(mean, sigma), obs=y)
        return mean


class BayesianRegressionModule(PyroBaseModuleClass):
    def __init__(self, **kwargs):

        super().__init__()
        self._model = BayesianRegressionPyroModel(**kwargs)
        self._guide = AutoNormal(
            self.model, init_loc_fn=init_to_mean, create_plates=self.model.create_plates
        )
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide

    @property
    def list_obs_plate_vars(self):
        return self.model.list_obs_plate_vars()


class BayesianRegressionModel(PyroSviTrainMixin, PyroSampleMixin, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        per_cell_weight=False,
    ):
        # in case any other model was created before that shares the same parameter names.
        clear_param_store()

        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        register_tensor_from_anndata(
            adata,
            registry_key="ind_x",
            adata_attr_name="obs",
            adata_key_name="_indices",
        )

        super().__init__(adata)

        self.module = BayesianRegressionModule(
            in_features=adata.shape[1],
            out_features=1,
            per_cell_weight=per_cell_weight,
        )
        self._model_summary_string = "BayesianRegressionModel"
        self.init_params_ = self._get_init_params(locals())


def test_pyro_bayesian_regression(save_path):
    use_gpu = int(torch.cuda.is_available())
    adata = synthetic_iid()
    # add index for each cell (provided to pyro plate for correct minibatching)
    adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
    register_tensor_from_anndata(
        adata,
        registry_key="ind_x",
        adata_attr_name="obs",
        adata_key_name="_indices",
    )
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(in_features=adata.shape[1], out_features=1)
    plan = PyroTrainingPlan(model)
    plan.n_obs_training = len(train_dl.indices)
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
    new_model = BayesianRegressionModule(in_features=adata.shape[1], out_features=1)
    # run model one step to get autoguide params
    try:
        new_model.load_state_dict(torch.load(model_save_path))
    except RuntimeError as err:
        if isinstance(new_model, PyroBaseModuleClass):
            plan = PyroTrainingPlan(new_model)
            plan.n_obs_training = len(train_dl.indices)
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
    # add index for each cell (provided to pyro plate for correct minibatching)
    adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
    register_tensor_from_anndata(
        adata,
        registry_key="ind_x",
        adata_attr_name="obs",
        adata_key_name="_indices",
    )
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(in_features=adata.shape[1], out_features=1)
    train_dl = AnnDataLoader(adata, shuffle=True, batch_size=128)
    plan = PyroTrainingPlan(model, loss_fn=pyro.infer.JitTrace_ELBO())
    plan.n_obs_training = len(train_dl.indices)
    trainer = Trainer(
        gpus=use_gpu, max_epochs=2, callbacks=[PyroJitGuideWarmup(train_dl)]
    )
    trainer.fit(plan, train_dl)

    # 100 features
    assert list(model.guide.state_dict()["locs.linear.weight_unconstrained"].shape) == [
        1,
        100,
    ]
    # 1 bias
    assert list(model.guide.state_dict()["locs.linear.bias_unconstrained"].shape) == [
        1,
    ]

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


def test_pyro_bayesian_train_sample_mixin():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    mod = BayesianRegressionModel(adata)
    mod.train(
        max_epochs=2,
        batch_size=128,
        plan_kwargs={"optim_kwargs": {"lr": 0.01}},
        use_gpu=use_gpu,
    )

    # 100 features
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=128, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10

    plt.plot(
        mod.history_["elbo_train"].index,
        mod.history_["elbo_train"].values.flatten(),
        label="elbo_train",
    )


def test_pyro_bayesian_train_sample_mixin_full_data():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    mod = BayesianRegressionModel(adata)
    mod.train(
        max_epochs=2,
        batch_size=None,
        plan_kwargs={"optim_kwargs": {"lr": 0.01}},
        use_gpu=use_gpu,
    )

    # 100 features
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=None, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10


def test_pyro_bayesian_train_sample_mixin_with_local():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    mod = BayesianRegressionModel(adata, per_cell_weight=True)
    mod.train(
        max_epochs=2,
        batch_size=128,
        plan_kwargs={"optim_kwargs": {"lr": 0.01}},
        train_size=1,  # does not work when there is a validation set.
        use_gpu=use_gpu,
    )

    # 100
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=128, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10
    assert samples["posterior_samples"]["per_cell_weights"].shape == (
        10,
        adata.n_obs,
        1,
    )


def test_pyro_bayesian_train_sample_mixin_with_local_full_data():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    mod = BayesianRegressionModel(adata, per_cell_weight=True)
    mod.train(
        max_epochs=2,
        batch_size=None,
        plan_kwargs={"optim_kwargs": {"lr": 0.01}},
        train_size=1,  # does not work when there is a validation set.
        use_gpu=use_gpu,
    )

    # 100
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=None, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10
    assert samples["posterior_samples"]["per_cell_weights"].shape == (
        10,
        adata.n_obs,
        1,
    )
