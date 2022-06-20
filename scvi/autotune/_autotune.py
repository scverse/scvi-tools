import inspect
import logging
import os
import warnings
from typing import Optional, Union, List

import anndata
import numpy as np
import ray
import scanpy as sc
import torch
from ray import tune
from ray.tune import CLIReporter, ExperimentAnalysis
from ray.tune.schedulers import ASHAScheduler
from ray.tune.schedulers.trial_scheduler import TrialScheduler

from scvi.model.base import BaseModelClass

from ._utils import apply_model_config, format_config, train_model

logger = logging.getLogger(__name__)


class Autotune:
    """
    Hyperparameter tuning using Ray Tune.
    Parameters
    ----------
    adata
        AnnData object we will tune the model on.
    model
        Model from scvi.model we will tune.
    num_epochs
        Number of epochs to tune the model over
    training_metrics
        Metrics to track during training.
    metric_functions
        For metrics calculated after training a model, like silhouette distance.
    model_hyperparams
        Config for the model hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    trainer_hyperparams
        Config for the trainer hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    plan_hyperparams
        Config for the training_plan hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    setup_args
        Non tunable anndata setup parameters
    test_effect_covariates
        Whether to test the effect of including covariates in the model. Default is `True`
    continuous_covariates
        continuous covariates from `adata` to tune the model on
    categorical_covariates
        Categorical covariates from `adata` to tune the model on
    test_effect_hvg
        Whether to test the effect of filtering adata with highly variable genes on model performance or not. Default is `False`.
    top_hvg
        Top hvgs to test. If `test_effect_hvg` is set, defaults to top 2000, 2500, 3000 and 3500 highly variable genes
    batch_key_hvg
        Column in adata.obs specifying the different batches in the data used for highly_variable_gene selection

     Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> tuner = Autotune(adata, scvi.model.SCVI, num_epochs=5)
    >>> best_model, analysis = tuner.run(metric="elbo_validation")
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        model: BaseModelClass,
        num_epochs: int = 2,
        training_metrics: Optional[List[str]] = None,
        metric_functions: Optional[dict] = None,
        model_hyperparams: Optional[dict] = None,
        trainer_hyperparams: Optional[dict] = None,
        plan_hyperparams: Optional[dict] = None,
        setup_args: Optional[dict] = None,
        continuous_covariates: Optional[List[str]] = None,
        categorical_covariates: Optional[List[str]] = None,
        test_effect_covariates: bool = False,
        test_effect_hvg: bool = False,
        top_hvg: Optional[List[int]] = None,
        batch_key_hvg: Optional[str] = None,
    ):
        training_metrics = training_metrics or []
        metric_functions = metric_functions or {}
        model_hyperparams = model_hyperparams or {}
        trainer_hyperparams = trainer_hyperparams or {}
        plan_hyperparams = plan_hyperparams or {}
        setup_args = setup_args or {}
        continuous_covariates = continuous_covariates or []
        categorical_covariates = categorical_covariates or []
        top_hvg = top_hvg or [None]

        if test_effect_hvg and not top_hvg:
            raise ValueError(
                "test_effect_hvg is set to True but no list of top_hvgs was given. Please"
                " provide a list"
            )

        self.adata = adata
        self.model = model
        self.setup_signature = inspect.signature(
            self.model.setup_anndata
        ).parameters  # Fetch available parameters
        self.training_metrics = training_metrics
        self.metric_functions = metric_functions
        self.model_hyperparams = model_hyperparams
        self.trainer_hyperparams = trainer_hyperparams
        self.plan_hyperparams = plan_hyperparams
        self.continious_covariates = {
            "continious_covariates": tune.choice(continious_covariates)
            if continious_covariates
        self.setup_args = setup_args

        continuous_covariates = [[covariate] for covariate in continuous_covariates]
        categorical_covariates = [[covariate] for covariate in categorical_covariates]
        if test_effect_covariates:
            continuous_covariates.append(None)
            categorical_covariates.append(None)

        self.continuous_covariates = {
            "continuous_covariate_keys": tune.choice(continuous_covariates)
            if continuous_covariates
            and "continuous_covariate_keys"
            in self.setup_signature  # Change default None only if setup_anndata signature has covariates
            else None
        }
        self.categorical_covariates = {
            "categorical_covariate_keys": tune.choice(categorical_covariates)
            if categorical_covariates
            and "categorical_covariate_keys"
            in self.setup_signature  # Change default None only if setup_anndata signature has covariates
            else None
        }

        if test_effect_hvg:
            self.test_effect_hvg = {"subset_adata": tune.choice([True, False])}
        else:
            self.test_effect_hvg = {"subset_adata": False}

        # Using conditional search spaces is not supported in all space search algorithms
        # If necessary we would need to re-factor to Optuna / HyperOptSearch spaces
        self.top_hvg = {
            "top_n": tune.sample_from(
                lambda spec: np.random.choice(top_hvg)
                if spec.config["subset_adata"]
                else None
            )
        }

        self.metrics = training_metrics
        self.reporter = CLIReporter(
            metric_columns=training_metrics + list(self.metric_functions.keys())
        )
        self.config = {}
        for att in [
            self.model_hyperparams,
            self.trainer_hyperparams,
            self.plan_hyperparams,
            self.continuous_covariates,
            self.categorical_covariates,
            self.test_effect_hvg,
            self.top_hvg,
        ]:
            if att is not None:
                self.config.update(att)
        self.batch_key = batch_key
        self.num_epochs = num_epochs


    def _trainable(self, config):

        model_config, trainer_config, plan_config, hvg_config = format_config(
            self, config
        )

        if "subset_adata" in hvg_config is True:
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat_v3",
                n_top_genes=hvg_config["top_n"],
                batch_key=self.batch_key_hvg,
                subset=True,
            )

        model = apply_model_config(self, model_config)
        train_model(self, model, trainer_config, plan_config)

    def run(
        self,
        metric: str = None,
        scheduler: TrialScheduler = None,
        mode: str = "min",
        name: str = "scvi-experiment",
        num_samples: int = 1,
        resources_per_trial: Optional[dict] = None,
        local_dir: str = "./ray_results",
        train: bool = False,
        **kwargs,
    ) -> Union[BaseModelClass, ExperimentAnalysis]:
        """
        Wrapper for `tune.run`. Searches for the configuration of model, trainer, and training_plan
        hyperparameters that minimize or maximize the provided metric.
        Parameters
        ----------
        metric
            Metric to optimize over in self.metrics or from self.training_funcs
        scheduler
            Ray tune scheduler for trials. If `None`, defaults to
            :class:`~ray.tune.schedulers.ASHAScheduler`, with `max_t` set to the number
            of epochs.
        mode
            "min" or "max" to maximize or minimize the objective metric
        name
            Name of this experiment.
        num_samples
            Number of times to sample hyperparameters from the configuration space
        resources_per_trial
            Dictionary specifying the number of `gpu` and `cpu` in the optimization
        local_dir
            Local dir to save training results to.
        train
            Whether to train the resulting best model. Defaults to `False`
        **kwargs
            Additional arguments for func:`ray.tune.run`

        Returns
        -------
        A tuple with the best model object and tune Analysis object
        """
        if not scheduler:
            scheduler = ASHAScheduler(
                max_t=self.num_epochs, grace_period=1, reduction_factor=2
            )

        if not resources_per_trial:
            if torch.cuda.is_available():
                resources_per_trial = {"cpu": 1, "gpu": 1}
            else:
                resources_per_trial = {"cpu": 1}

        if resources_per_trial["cpu"] > 1 and resources_per_trial["gpu"]:
            resources_per_trial["gpu"] /= resources_per_trial["cpu"]

        if (
            metric in self.metrics
            and (
                "continuous_covariate_keys"
                or "categorical_covariate_keys" in self.setup_args
            )
            or self.top_hvg
        ):
            warnings.warn(
                "You are optimizing over {} and testing different model architectures. "
                "This is not a recommended approach as the metric will be influenced"
                " by the architecture and not the model hyperparameters. Consider "
                "optimizing over a non-training metric, such as `autotune.metrics.silhouette_score`".format(
                    metric
                )
            )

        ray.init(num_cpus=resources_per_trial["cpu"], ignore_reinit_error=True)
        analysis = tune.run(
            self._trainable,
            metric=metric,
            mode=mode,
            config=self.config,
            num_samples=num_samples,
            scheduler=scheduler,
            progress_reporter=self.reporter,
            name=name,
            resources_per_trial=resources_per_trial,
            local_dir=local_dir,
            raise_on_failed_trial=False,
            **kwargs,
        )
        ray.shutdown()
        logger.info(
            "Best hyperparameters found were: {}".format(analysis.get_best_config())
        )

        best_config = analysis.best_config

        _, trainer_config, plan_config, hvg_config = format_config(self, best_config)

        if "subset_adata" in hvg_config is True:
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat_v3",
                n_top_genes=hvg_config["top_n"],
                batch_key=self.batch_key_hvg,
                subset=True,
            )

        # retrieve and load best model
        best_checkpoint = analysis.best_checkpoint
        best_model = self.model.load(
            dir_path=os.path.join(best_checkpoint, "checkpoint"), adata=self.adata
        )
        best_model = train_model(
            self,
            best_model,
            trainer_config,
            plan_config,
            return_model=True,
            callbacks=False,
        )

        return best_model, analysis
