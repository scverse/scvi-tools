from typing import Tuple, Union

from scvi.model.base import BaseModelClass

from ._callbacks import ModelSave, _TuneReportMetricFunctionsCallback


def format_config(self, config: dict = None) -> Tuple[dict, ...]:
    """
    Format tune's config dictionaries to use as input in scvi-tools' workflow

    ----------
    config
        Dictionary specifying the search space used by `ray.tune.run`

    Returns
    -------
    A tuple of dictionaries specifying the hyperparameter configuration
    """
    model_config = {}
    trainer_config = {}
    plan_config = {}
    hvg_config = {}
    for key in config:
        if key in self.model_hyperparams:
            model_config[key] = config[key]
        elif key in self.trainer_hyperparams:
            trainer_config[key] = config[key]
        elif key in self.plan_hyperparams:
            plan_config[key] = config[key]
        elif key in self.continuous_covariates:
            self.setup_args[key] = config[
                key
            ]  # store tunable covariates in setup_args dict so we can generalize to other models
        elif key in self.categorical_covariates:
            self.setup_args[key] = config[key]
        elif key in self.test_effect_hvg:
            hvg_config[key] = config[key]
        elif key in self.top_hvg:
            hvg_config[key] = config[key]

    return model_config, trainer_config, plan_config, hvg_config


def apply_model_config(
    self,
    model_config: dict = None,
) -> BaseModelClass:
    """
    Initializes the specified model in Autotune.model with a given anndata setup and model and training hyperparameters.

    ----------
    model_config
        Dictionary specifying the model hyperparameters used to initialize the model

    Returns
    -------
    An object of class :class:`~scvi.model.base.BaseModelClass`
    """

    self.model.setup_anndata(
        self.adata,
        **self.setup_args,
    )
    model = self.model(self.adata, **model_config)
    return model


def train_model(
    self,
    model: BaseModelClass,
    trainer_config: dict = None,
    plan_config: dict = None,
    return_model: bool = False,
    callbacks: bool = True,
) -> Union[None, BaseModelClass]:
    """
    Trains the specified model in Autotune.model with a given training plan.

    ----------
    model_config
        Dictionary specifying the model hyperparameters used to initialize the model
    trainer_config
        Dictionary specifying the trainer hyperparameters passed to :class::`~scvi.train.Trainer`
    plan_config
        Dictionary specifying the plan configuration passed to :func:`~scvi.train.TrainingPlan`
    return_model
        Boolean specifying whether to return the trained model or not

    Returns
    -------
    If return_model is set to True, a trained object of class :class:`~scvi.model.base.BaseModelClass`
    """
    if not callbacks:
        model.train(
            **trainer_config, plan_kwargs=plan_config, max_epochs=self.num_epochs
        )
    else:
        model.train(
            **trainer_config,
            plan_kwargs=plan_config,
            callbacks=[
                ModelSave(model),
                _TuneReportMetricFunctionsCallback(
                    metrics=self.metrics,
                    on="validation_end",
                    model=model,
                    metric_functions=self.metric_functions,
                ),
            ],
            check_val_every_n_epoch=1,
            max_epochs=self.num_epochs,
        )
    return model if return_model else None
