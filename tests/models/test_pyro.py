import os
from typing import Optional

import numpy as np
import pyro
import pyro.distributions as dist
import torch
from anndata import AnnData
from pyro import clear_param_store
from pyro.infer.autoguide import AutoNormal, init_to_mean
from pyro.nn import PyroModule, PyroSample
from torch import nn

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager, synthetic_iid
from scvi.data.fields import CategoricalObsField, LayerField, NumericalObsField
from scvi.dataloaders import AnnDataLoader
from scvi.model import AmortizedLDA
from scvi.model.base import (
    BaseModelClass,
    PyroJitGuideWarmup,
    PyroModelGuideWarmup,
    PyroSampleMixin,
    PyroSviTrainMixin,
)
from scvi.module.base import PyroBaseModuleClass
from scvi.nn import DecoderSCVI, Encoder
from scvi.train import LowLevelPyroTrainingPlan, PyroTrainingPlan, Trainer


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
        # Using the lambda here means that Pyro recreates the prior every time
        # it retrieves it. In this case this is first when the model is called, which is
        # also when the auto guide params are created.
        # This effectively allows these priors to move with the model's device with the
        # expense of dynamic recreation.
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
        """
        Function for creating plates is needed when using AutoGuides.

        Should have the same call signature as model.
        """
        return pyro.plate("obs_plate", size=self.n_obs, dim=-2, subsample=ind_x)

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
        return {
            "name": "obs_plate",
            "in": [0],  # model args index for expression data
            "sites": {"per_cell_weights": 1},
        }

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        y = tensor_dict[REGISTRY_KEYS.LABELS_KEY]
        ind_x = tensor_dict[REGISTRY_KEYS.INDICES_KEY].long().squeeze()
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

        super().__init__(adata)

        self.module = BayesianRegressionModule(
            in_features=adata.shape[1],
            out_features=1,
            per_cell_weight=per_cell_weight,
        )
        self._model_summary_string = "BayesianRegressionModel"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> Optional[AnnData]:
        setup_method_args = cls._get_setup_method_args(**locals())

        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, None, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, None),
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


def _create_indices_adata_manager(adata: AnnData) -> AnnDataManager:
    # add index for each cell (provided to pyro plate for correct minibatching)
    adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
    anndata_fields = [
        LayerField(REGISTRY_KEYS.X_KEY, None, is_count_data=True),
        CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, None),
        NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
    ]
    adata_manager = AnnDataManager(fields=anndata_fields)
    adata_manager.register_fields(adata)
    return adata_manager


def test_pyro_bayesian_regression_low_level():
    use_gpu = int(torch.cuda.is_available())
    adata = synthetic_iid()
    adata_manager = _create_indices_adata_manager(adata)
    train_dl = AnnDataLoader(adata_manager, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(in_features=adata.shape[1], out_features=1)
    plan = LowLevelPyroTrainingPlan(model)
    plan.n_obs_training = len(train_dl.indices)
    trainer = Trainer(
        gpus=use_gpu,
        max_epochs=2,
        callbacks=[PyroModelGuideWarmup(train_dl)],
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


def test_pyro_bayesian_regression(save_path):
    use_gpu = int(torch.cuda.is_available())
    adata = synthetic_iid()
    adata_manager = _create_indices_adata_manager(adata)
    train_dl = AnnDataLoader(adata_manager, shuffle=True, batch_size=128)
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
    adata_manager = _create_indices_adata_manager(adata)
    train_dl = AnnDataLoader(adata_manager, shuffle=True, batch_size=128)
    pyro.clear_param_store()
    model = BayesianRegressionModule(in_features=adata.shape[1], out_features=1)
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


def test_pyro_bayesian_save_load(save_path):
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    BayesianRegressionModel.setup_anndata(adata)
    mod = BayesianRegressionModel(adata)
    mod.train(
        max_epochs=2,
        batch_size=128,
        lr=0.01,
        use_gpu=use_gpu,
    )

    mod.module.cpu()
    quants = mod.module.guide.quantiles([0.5])
    sigma_median = quants["sigma"][0].detach().cpu().numpy()
    linear_median = quants["linear.weight"][0].detach().cpu().numpy()

    model_save_path = os.path.join(save_path, "test_pyro_bayesian/")
    mod.save(model_save_path)

    # Test setting `on_load_kwargs`
    mod.module.on_load_kwargs = {"batch_size": 8}
    mod = BayesianRegressionModel.load(model_save_path, adata=adata)

    quants = mod.module.guide.quantiles([0.5])
    sigma_median_new = quants["sigma"][0].detach().cpu().numpy()
    linear_median_new = quants["linear.weight"][0].detach().cpu().numpy()

    np.testing.assert_array_equal(sigma_median_new, sigma_median)
    np.testing.assert_array_equal(linear_median_new, linear_median)


def test_pyro_bayesian_train_sample_mixin():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    BayesianRegressionModel.setup_anndata(adata)
    mod = BayesianRegressionModel(adata)
    mod.train(
        max_epochs=2,
        batch_size=128,
        lr=0.01,
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


def test_pyro_bayesian_train_sample_mixin_full_data():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    BayesianRegressionModel.setup_anndata(adata)
    mod = BayesianRegressionModel(adata)
    mod.train(
        max_epochs=2,
        batch_size=None,
        lr=0.01,
        use_gpu=use_gpu,
    )

    # 100 features
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=adata.n_obs, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10


def test_pyro_bayesian_train_sample_mixin_with_local():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    BayesianRegressionModel.setup_anndata(adata)
    mod = BayesianRegressionModel(adata, per_cell_weight=True)
    mod.train(
        max_epochs=2,
        batch_size=128,
        lr=0.01,
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


def test_pyro_bayesian_train_sample_mixin_with_local_full_data():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    BayesianRegressionModel.setup_anndata(adata)
    mod = BayesianRegressionModel(adata, per_cell_weight=True)
    mod.train(
        max_epochs=2,
        batch_size=None,
        lr=0.01,
        train_size=1,  # does not work when there is a validation set.
        use_gpu=use_gpu,
    )

    # 100
    assert list(
        mod.module.guide.state_dict()["locs.linear.weight_unconstrained"].shape
    ) == [1, 100]

    # test posterior sampling
    samples = mod.sample_posterior(
        num_samples=10, use_gpu=use_gpu, batch_size=adata.n_obs, return_samples=True
    )

    assert len(samples["posterior_samples"]["sigma"]) == 10
    assert samples["posterior_samples"]["per_cell_weights"].shape == (
        10,
        adata.n_obs,
        1,
    )


class FunctionBasedPyroModule(PyroBaseModuleClass):
    def __init__(self, n_input: int, n_latent: int, n_hidden: int, n_layers: int):

        super().__init__()
        self.n_input = n_input
        self.n_latent = n_latent
        self.epsilon = 5.0e-3
        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        self.encoder = Encoder(
            n_input,
            n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0.1,
            return_dist=True,
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        self.decoder = DecoderSCVI(
            n_latent,
            n_input,
            n_layers=n_layers,
            n_hidden=n_hidden,
        )
        # This gene-level parameter modulates the variance of the observation distribution
        self.px_r = torch.nn.Parameter(torch.ones(self.n_input))

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        log_library = torch.log(torch.sum(x, dim=1, keepdim=True) + 1e-6)
        return (x, log_library), {}

    def model(self, x, log_library):
        # register PyTorch module `decoder` with Pyro
        pyro.module("scvi", self)
        with pyro.plate("data", x.shape[0]):
            # setup hyperparameters for prior p(z)
            z_loc = x.new_zeros(torch.Size((x.shape[0], self.n_latent)))
            z_scale = x.new_ones(torch.Size((x.shape[0], self.n_latent)))
            # sample from prior (value will be sampled by guide when computing the ELBO)
            z = pyro.sample("latent", dist.Normal(z_loc, z_scale).to_event(1))
            # decode the latent code z
            px_scale, _, px_rate, px_dropout = self.decoder("gene", z, log_library)
            # build count distribution
            nb_logits = (px_rate + self.epsilon).log() - (
                self.px_r.exp() + self.epsilon
            ).log()
            x_dist = dist.ZeroInflatedNegativeBinomial(
                gate_logits=px_dropout, total_count=self.px_r.exp(), logits=nb_logits
            )
            # score against actual counts
            pyro.sample("obs", x_dist.to_event(1), obs=x)

    def guide(self, x, log_library):
        # define the guide (i.e. variational distribution) q(z|x)
        pyro.module("scvi", self)
        with pyro.plate("data", x.shape[0]):
            # use the encoder to get the parameters used to define q(z|x)
            x_ = torch.log(1 + x)
            qz, _ = self.encoder(x_)
            # sample the latent code z
            pyro.sample("latent", dist.Normal(qz.loc, qz.scale).to_event(1))


class FunctionBasedPyroModel(PyroSviTrainMixin, PyroSampleMixin, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
    ):
        # in case any other model was created before that shares the same parameter names.
        clear_param_store()

        super().__init__(adata)

        self.module = FunctionBasedPyroModule(
            n_input=adata.n_vars,
            n_hidden=32,
            n_latent=5,
            n_layers=1,
        )
        self._model_summary_string = "FunctionBasedPyroModel"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> Optional[AnnData]:
        setup_method_args = cls._get_setup_method_args(**locals())

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, None, is_count_data=True),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


def test_function_based_pyro_module():
    use_gpu = torch.cuda.is_available()
    adata = synthetic_iid()
    FunctionBasedPyroModel.setup_anndata(adata)
    mod = FunctionBasedPyroModel(adata)
    mod.train(
        max_epochs=1,
        batch_size=256,
        lr=0.01,
        use_gpu=use_gpu,
    )


def test_lda_model_single_step():
    n_topics = 5
    adata = synthetic_iid()
    AmortizedLDA.setup_anndata(adata)
    mod1 = AmortizedLDA(
        adata, n_topics=n_topics, cell_topic_prior=1.5, topic_feature_prior=1.5
    )
    mod1.train(max_steps=1, max_epochs=10)
    assert len(mod1.history["elbo_train"]) == 1


def test_lda_model():
    use_gpu = torch.cuda.is_available()
    n_topics = 5
    adata = synthetic_iid()

    # Test with float and Sequence priors.
    AmortizedLDA.setup_anndata(adata)
    mod1 = AmortizedLDA(
        adata, n_topics=n_topics, cell_topic_prior=1.5, topic_feature_prior=1.5
    )
    mod1.train(
        max_epochs=1,
        batch_size=256,
        lr=0.01,
        use_gpu=use_gpu,
    )
    mod2 = AmortizedLDA(
        adata,
        n_topics=n_topics,
        cell_topic_prior=[1.5 for _ in range(n_topics)],
        topic_feature_prior=[1.5 for _ in range(adata.n_vars)],
    )
    mod2.train(
        max_epochs=1,
        batch_size=256,
        lr=0.01,
        use_gpu=use_gpu,
    )

    mod = AmortizedLDA(adata, n_topics=n_topics)
    mod.train(
        max_epochs=5,
        batch_size=256,
        lr=0.01,
        use_gpu=use_gpu,
    )
    adata_gbt = mod.get_feature_by_topic().to_numpy()
    assert np.allclose(adata_gbt.sum(axis=0), 1)
    adata_lda = mod.get_latent_representation(adata).to_numpy()
    assert (
        adata_lda.shape == (adata.n_obs, n_topics)
        and np.all((adata_lda <= 1) & (adata_lda >= 0))
        and np.allclose(adata_lda.sum(axis=1), 1)
    )
    mod.get_elbo()
    mod.get_perplexity()

    adata2 = synthetic_iid()
    AmortizedLDA.setup_anndata(adata2)
    adata2_lda = mod.get_latent_representation(adata2).to_numpy()
    assert (
        adata2_lda.shape == (adata2.n_obs, n_topics)
        and np.all((adata2_lda <= 1) & (adata2_lda >= 0))
        and np.allclose(adata2_lda.sum(axis=1), 1)
    )
    mod.get_elbo(adata2)
    mod.get_perplexity(adata2)


def test_lda_model_save_load(save_path):
    use_gpu = torch.cuda.is_available()
    n_topics = 5
    adata = synthetic_iid()
    AmortizedLDA.setup_anndata(adata)
    mod = AmortizedLDA(adata, n_topics=n_topics)
    mod.train(
        max_epochs=5,
        batch_size=256,
        lr=0.01,
        use_gpu=use_gpu,
    )
    hist_elbo = mod.history_["elbo_train"]

    feature_by_topic_1 = mod.get_feature_by_topic(n_samples=5000)
    latent_1 = mod.get_latent_representation(n_samples=6000)

    save_path = os.path.join(save_path, "tmp")
    mod.save(save_path, overwrite=True, save_anndata=True)
    mod = AmortizedLDA.load(save_path)

    np.testing.assert_array_equal(mod.history_["elbo_train"], hist_elbo)

    feature_by_topic_2 = mod.get_feature_by_topic(n_samples=5000)
    latent_2 = mod.get_latent_representation(n_samples=6000)
    np.testing.assert_almost_equal(
        feature_by_topic_1.to_numpy(), feature_by_topic_2.to_numpy(), decimal=2
    )
    np.testing.assert_almost_equal(latent_1.to_numpy(), latent_2.to_numpy(), decimal=2)
