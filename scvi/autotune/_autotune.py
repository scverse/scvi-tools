import logging
from typing import Optional, Union

import anndata
<<<<<<< HEAD
=======
import ray
import scanpy as sc
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
import torch
from ray import tune
from ray.tune import CLIReporter, ExperimentAnalysis
from ray.tune.schedulers import ASHAScheduler
from ray.tune.schedulers.trial_scheduler import TrialScheduler
<<<<<<< HEAD
from ray.tune.suggest.hyperopt import HyperOptSearch
from ray.tune.suggest.search import SearchAlgorithm
=======
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)

from scvi.model.base import BaseModelClass

from ._callbacks import ModelSave, _TuneReportMetricFunctionsCallback

logger = logging.getLogger(__name__)


class Autotune:
    """
    Hyperparameter tuning using Ray Tune.
<<<<<<< HEAD

=======
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    Parameters
    ----------
    adata
        AnnData object we will tune the model on.
    model
        Model from scvi.model we will tune.
<<<<<<< HEAD
=======
    num_epochs
        Number of epochs to tune the model over
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
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
<<<<<<< HEAD
=======
    test_effect_covariates
        Whether to test the effect of including covariates in the model. Default is `True`
    continious_covariates
        Continious covariates from `adata` to tune the model on
    categorical_covariates
        Categorical covariates from `adata` to tune the model on
    test_effect_hvg
        Whether to test the effect of filtering adata with highly variable genes on model performance or not. Default is `False`.
    top_hvg
        Top hvgs to test. If `test_effect_hvg` is set, defaults to top 2000, 2500, 3000 and 3500 highly variable genes
    batch_key
        Column in adata.obs specifying the different batches in the data.

     Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> tuner = Autotune(adata, scvi.model.SCVI, num_epochs=5)
    >>> best_model, analysis = tuner.run(metric="elbo_validation")
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        model: BaseModelClass,
<<<<<<< HEAD
        training_metrics: Optional[list] = None,
=======
        num_epochs: int = 2,
        training_metrics: Optional[list[str]] = None,
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        metric_functions: Optional[dict] = None,
        model_hyperparams: Optional[dict] = None,
        trainer_hyperparams: Optional[dict] = None,
        plan_hyperparams: Optional[dict] = None,
<<<<<<< HEAD
        num_epochs: int = 2,
=======
        continious_covariates: Optional[list[str]] = None,
        categorical_covariates: Optional[list[str]] = None,
        test_effect_covariates: bool = False,
        batch_key: Optional[str] = None,
        test_effect_hvg: bool = False,
        top_hvg: Optional[list[int]] = None,
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
    ):
        if not training_metrics:
            training_metrics = []
        if not metric_functions:
            metric_functions = {}
        if not model_hyperparams:
            model_hyperparams = {}
        if not trainer_hyperparams:
            trainer_hyperparams = {}
        if not plan_hyperparams:
            plan_hyperparams = {}
<<<<<<< HEAD
=======
        if not continious_covariates:
            continious_covariates = []
        if not categorical_covariates:
            categorical_covariates = []
        if not top_hvg:
            top_hvg = []
        if test_effect_hvg and not top_hvg:
            raise ValueError(
                "test_effect_hvg is set to True but no list of top_hvgs was given. Please"
                " provide a list"
            )

        continious_covariates = [[covariate] for covariate in continious_covariates]
        categorical_covariates = [[covariate] for covariate in categorical_covariates]
        if test_effect_covariates:
            continious_covariates.append(None)
            categorical_covariates.append(None)

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        self.adata = adata
        self.model = model
        self.training_metrics = training_metrics
        self.metric_functions = metric_functions
        self.model_hyperparams = model_hyperparams
        self.trainer_hyperparams = trainer_hyperparams
        self.plan_hyperparams = plan_hyperparams
<<<<<<< HEAD
=======
        self.continious_covariates = {
            "continious_covariates": tune.choice(continious_covariates)
            if continious_covariates
            else None
        }
        self.categorical_covariates = {
            "categorical_covariates": tune.choice(categorical_covariates)
            if categorical_covariates
            else None
        }

        if test_effect_hvg:
            self.test_effect_hvg = {"subset_adata": tune.choice([True, False])}
        else:
            self.test_effect_hvg = {"subset_adata": False}

        self.top_hvg = {
            "top_n": tune.choice(top_hvg) if test_effect_hvg or top_hvg else None
        }

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        self.metrics = training_metrics
        self.reporter = CLIReporter(
            metric_columns=training_metrics + list(self.metric_functions.keys())
        )
        self.config = {}
<<<<<<< HEAD
        for d in [model_hyperparams, trainer_hyperparams, plan_hyperparams]:
            if d is not None:
                self.config.update(d)
=======
        for att in [
            self.model_hyperparams,
            self.trainer_hyperparams,
            self.plan_hyperparams,
            self.continious_covariates,
            self.categorical_covariates,
            self.test_effect_hvg,
            self.top_hvg,
        ]:
            if att is not None:
                self.config.update(att)
        self.batch_key = batch_key
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        self.num_epochs = num_epochs

    def _trainable(self, config):
        model_config = {}
        trainer_config = {}
        plan_config = {}
<<<<<<< HEAD
=======
        continious_covariate_config = {}
        categorical_covariate_config = {}
        hvg_config = {}

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        for key in config:
            if key in self.model_hyperparams:
                model_config[key] = config[key]
            elif key in self.trainer_hyperparams:
                trainer_config[key] = config[key]
            elif key in self.plan_hyperparams:
                plan_config[key] = config[key]
<<<<<<< HEAD

=======
            elif key in self.continious_covariates:
                continious_covariate_config[key] = config[key]
            elif key in self.categorical_covariates:
                categorical_covariate_config[key] = config[key]
            elif (
                (self.test_effect_hvg and self.top_hvg) is not None
                and key in self.test_effect_hvg
                or self.top_hvg
            ):
                hvg_config[key] = config[key]

        if "subset_adata" in hvg_config is True:
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat_v3",
                n_top_genes=hvg_config["top_n"],
                batch_key=self.batch_key,
                subset=True,
            )
        self.model.setup_anndata(
            self.adata,
            batch_key=self.batch_key,
            continuous_covariate_keys=continious_covariate_config[
                "continious_covariates"
            ],
            categorical_covariate_keys=categorical_covariate_config[
                "categorical_covariates"
            ],
        )
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        _model = self.model(self.adata, **model_config)
        _model.train(
            **trainer_config,
            plan_kwargs=plan_config,
            callbacks=[
                ModelSave(_model),
                _TuneReportMetricFunctionsCallback(
                    metrics=self.metrics,
                    on="validation_end",
                    model=_model,
                    metric_functions=self.metric_functions,
                ),
            ],
            check_val_every_n_epoch=1,
            max_epochs=self.num_epochs,
        )

    def run(
        self,
        metric: str = None,
        scheduler: TrialScheduler = None,
<<<<<<< HEAD
        search_alg: SearchAlgorithm = None,
        mode: str = "min",
        name: str = "scvi-experiment",
        num_samples: int = 10,
        resources_per_trial: Optional[dict] = None,
        local_dir: str = "./ray_results",
=======
        mode: str = "min",
        name: str = "scvi-experiment",
        num_samples: int = 2,
        resources_per_trial: Optional[dict] = None,
        local_dir: str = "./ray_results",
        train: bool = False,
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        **kwargs,
    ) -> Union[BaseModelClass, ExperimentAnalysis]:
        """
        Wrapper for `tune.run`. Searches for the configuration of model, trainer, and training_plan
        hyperparameters that minimize or maximize the provided metric.
<<<<<<< HEAD

=======
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        Parameters
        ----------
        metric
            Metric to optimize over in self.metrics or from self.training_funcs
        scheduler
            Ray tune scheduler for trials. If `None`, defaults to
<<<<<<< HEAD
            :class:`~ray.tune.schedulers.ASHAScheduler`.
        search_alg
            Search algorithm from `tune.suggest`. If `None`, defaults to
            :class:`~ray.tune.suggest.hyperopt.HyperOptSearch`.
=======
            :class:`~ray.tune.schedulers.ASHAScheduler`, with `max_t` set to the number
            of epochs.
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        mode
            "min" or "max" to maximize or minimize the objective metric
        name
            Name of this experiment.
        num_samples
<<<<<<< HEAD
            Number of times to sample hyperparameters from the configuration space.
        local_dir
            Local dir to save training results to.
=======
            Number of times to sample hyperparameters from the configuration space
        resources_per_trial
            Dictionary specifying the number of `gpu` and `cpu` in the optimization
        local_dir
            Local dir to save training results to.
        train
            Whether to train the resulting best model. Defaults to `False`
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        **kwargs
            Additional arguments for func:`ray.tune.run`

        Returns
        -------
        A tuple with the best model object and tune Analysis object
<<<<<<< HEAD

        """
        if not scheduler:
            scheduler = ASHAScheduler(max_t=2, grace_period=1, reduction_factor=2)

        if not search_alg:
            # metric, mode, and config will be passed in tune.run
            search_alg = HyperOptSearch()
=======
        """
        if not scheduler:
            scheduler = ASHAScheduler(
                max_t=self.num_epochs, grace_period=1, reduction_factor=2
            )
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)

        if not resources_per_trial:
            if torch.cuda.is_available():
                resources_per_trial = {"cpu": 1, "gpu": 1}
            else:
                resources_per_trial = {"cpu": 1}

<<<<<<< HEAD
=======
        if resources_per_trial["cpu"] > 1 and resources_per_trial["gpu"]:
            resources_per_trial["gpu"] /= resources_per_trial["cpu"]

        ray.init(num_cpus=resources_per_trial["cpu"], ignore_reinit_error=True)
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        analysis = tune.run(
            self._trainable,
            metric=metric,
            mode=mode,
            config=self.config,
            num_samples=num_samples,
            scheduler=scheduler,
<<<<<<< HEAD
            search_alg=search_alg,
=======
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
            progress_reporter=self.reporter,
            name=name,
            resources_per_trial=resources_per_trial,
            local_dir=local_dir,
<<<<<<< HEAD
            **kwargs,
        )
        best_config = analysis.best_config
        logger.info("Best hyperparameters found were: ", best_config)
        # Get the checkpoint path of the best trial of the experiment
        model_config = {}
        trainer_config = {}
        plan_config = {}
=======
            raise_on_failed_trial=False,
            **kwargs,
        )
        ray.shutdown()
        logger.info(
            "Best hyperparameters found were: {}".format(analysis.get_best_config())
        )
        best_config = analysis.best_config
        model_config = {}
        trainer_config = {}
        plan_config = {}
        continious_covariate_config = {}
        categorical_covariate_config = {}
        hvg_config = {}
>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        for key in best_config:
            if key in self.model_hyperparams:
                model_config[key] = best_config[key]
            elif key in self.trainer_hyperparams:
                trainer_config[key] = best_config[key]
            elif key in self.plan_hyperparams:
                plan_config[key] = best_config[key]
<<<<<<< HEAD

        # retrieve and load best model
        best_checkpoint = analysis.best_checkpoint
        best_model = self.model(self.adata, **model_config)
        best_model.load(adata=self.adata, dir_path=best_checkpoint + "checkpoint")
=======
            elif key in self.continious_covariates:
                continious_covariate_config[key] = best_config[key]
            elif key in self.categorical_covariates:
                categorical_covariate_config[key] = best_config[key]
            elif (
                (self.test_effect_hvg and self.top_hvg) is not None
                and key in self.test_effect_hvg
                or self.top_hvg
            ):
                hvg_config[key] = best_config[key]

        # retrieve and load best model
        best_checkpoint = analysis.best_checkpoint

        if "subset_adata" in hvg_config is True:
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat_v3",
                n_top_genes=hvg_config["top_n"],
                batch_key=self.batch_key,
                subset=True,
            )

        self.model.setup_anndata(
            self.adata,
            batch_key=self.batch_key,
            continuous_covariate_keys=continious_covariate_config[
                "continious_covariates"
            ],
            categorical_covariate_keys=categorical_covariate_config[
                "categorical_covariates"
            ],
        )

        best_model = self.model(self.adata, **model_config)
        best_model.load(adata=self.adata, dir_path=best_checkpoint + "checkpoint")
        if train:
            best_model.train(
                **trainer_config,
                plan_kwargs=plan_config,
                check_val_every_n_epoch=1,
                max_epochs=self.num_epochs,
            )

>>>>>>> 535fc45f (Ray.tune for parameter optimization. Skeleton based on branch: michael/autotune. Included more funcitonality)
        return best_model, analysis
