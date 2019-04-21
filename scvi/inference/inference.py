import copy
from collections import defaultdict
from functools import partial
from typing import Any, Dict, Type, Union

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK
from hyperopt.mongoexp import MongoTrials

import matplotlib.pyplot as plt
import torch

from . import Trainer
from ..dataset import GeneExpressionDataset
from ..models import VAE

plt.switch_backend('agg')


class UnsupervisedTrainer(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\*\*kwargs: Other keywords arguments from the general Trainer class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, test_size=None, kl=None, **kwargs):
        super().__init__(model, gene_dataset, **kwargs)
        self.kl = kl
        if type(self) is UnsupervisedTrainer:
            self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size, test_size)
            self.train_set.to_monitor = ['ll']
            self.test_set.to_monitor = ['ll']

    @property
    def posteriors_loop(self):
        return ['train_set']

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / 400)  # self.n_epochs)


class AdapterTrainer(UnsupervisedTrainer):
    def __init__(self, model, gene_dataset, posterior_test, frequency=5):
        super().__init__(model, gene_dataset, frequency=frequency)
        self.test_set = posterior_test
        self.test_set.to_monitor = ['ll']
        self.params = list(self.model.z_encoder.parameters()) + list(self.model.l_encoder.parameters())
        self.z_encoder_state = copy.deepcopy(model.z_encoder.state_dict())
        self.l_encoder_state = copy.deepcopy(model.l_encoder.state_dict())

    @property
    def posteriors_loop(self):
        return ['test_set']

    def train(self, n_path=10, n_epochs=50, **kwargs):
        for i in range(n_path):
            # Re-initialize to create new path
            self.model.z_encoder.load_state_dict(self.z_encoder_state)
            self.model.l_encoder.load_state_dict(self.l_encoder_state)
            super().train(n_epochs, params=self.params, **kwargs)

        return min(self.history["ll_test_set"])


def auto_tuned_scvi_model(
    gene_dataset: GeneExpressionDataset,
    model_class: VAE = VAE,
    trainer_class: Trainer = UnsupervisedTrainer,
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    space: dict = None,
    use_batches: bool = False,
    max_evals: int = 10,
    parallel: bool = False,
    verbose: bool = True,
) -> (Type[Trainer], Trials):
    """Perform automatic hyperparameter optimization of an scVI model
    and return best model and hyperopt Trials object.
    Trials object contains hyperparameter space and loss history for each trial.
    We provide a default hyperparameter search space (see source code),
    but we recommend the user to build a custom one for each application.
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).

    :param gene_dataset: scVI gene dataset
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: Trainer class (e.g ``UnsupervisedTrainer``)
    :param model_specific_kwargs: dict of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: dict of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param space: dict containing up to three sub-dicts with keys 'model_tunable_kwargs',
    'trainer_tunable_kwargs' or 'train_func_tunable_kwargs'.
    Each of those dict contains hyperopt defined parameter spaces (e.g. ``hp.choice(..)``)
    which will be passed to the corresponding object : model, trainer or train method
    when performing hyperoptimization. Default: mutable, see source code.
    :param use_batches: If False, pass n_batch=0 to model else pass gene_dataset.n_batches
    :param max_evals: Maximum number of trainings.
    :return:

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tuned_scvi_parameters(gene_dataset)
        >>> posterior = best_trainer.create_posterior()
    """

    # default specific kwargs
    model_specific_kwargs = model_specific_kwargs if model_specific_kwargs else {}
    trainer_specific_kwargs = trainer_specific_kwargs if trainer_specific_kwargs else {}
    train_func_specific_kwargs = train_func_specific_kwargs if train_func_specific_kwargs else {}

    # default early stopping
    if 'early_stopping_kwargs' not in trainer_specific_kwargs:
        early_stopping_kwargs = {
            'early_stopping_metric': 'll',
            'save_best_state_metric': 'll',
            'patience': 15,
            'threshold': 3,
        }
        trainer_specific_kwargs['early_stopping_kwargs'] = early_stopping_kwargs

    # default search space
    # FIXME: what should the defaults be?
    space = space if space else {
        'model_tunable_kwargs': {
            'n_latent': 5 + hp.randint('n_latent', 5),
            'n_hidden': hp.choice('n_hidden', [64, 128, 256]),
            'n_layers': 1 + hp.randint('n_layers', 3),
            'dropout_rate': hp.uniform('dropout_rate', 0.1, 0.9),
        },
        'trainer_tunable_kwargs': {
            'kl': hp.uniform('kl', 0.5, 1.5),
        },
        'train_func_tunable_kwargs': {
            'n_epochs': hp.choice('n_epochs', [100, 200]),
            'lr': hp.choice('lr', [0.01, 0.001, 0.0001]),
        },
    }

    if verbose:
        print('Fixed parameters: ')
        print('model:')
        print(model_specific_kwargs)
        print('trainer: ')
        print(trainer_specific_kwargs)
        print('train method: ')
        print(train_func_specific_kwargs)

    # build a partial objective function restricted to the search space
    objective_hyperopt = partial(
        _objective_function,
        **{
            'gene_dataset': gene_dataset,
            'model_class': model_class,
            'trainer_class': trainer_class,
            'model_specific_kwargs': model_specific_kwargs,
            'trainer_specific_kwargs': trainer_specific_kwargs,
            'train_func_specific_kwargs': train_func_specific_kwargs,
            'use_batches': use_batches,
            'verbose': verbose,
        }
    )

    # run hyperoptimization
    trials = MongoTrials(
        'mongo://localhost:1234/scvi_db/jobs',
        exp_key='exp1'
    ) if parallel else Trials()

    _ = fmin(
        objective_hyperopt,
        space=space,
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials,
    )

    # return best model, trained
    best_space = trials.best_trial['result']['space']
    best_trainer = objective_hyperopt(best_space, is_best_training=True)

    return best_trainer, trials


def _objective_function(
    space: dict,
    gene_dataset: GeneExpressionDataset,
    model_class: Type[VAE] = VAE,
    trainer_class: Type[Trainer] = UnsupervisedTrainer,
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    use_batches: bool = False,
    verbose: bool = True,
    is_best_training: bool = False,
) -> Union[Dict[str, Any], Type[VAE]]:
    """Objective function for automatic hyperparameter optimization.
    Train a scVI model and return the best value of the early-stopping metric (e.g, log-likelihood).
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).

    :param space: dict containing up to three sub-dicts with keys 'model_tunable_kwargs',
    'trainer_tunable_kwargs' or 'train_func_tunable_kwargs'.
    Each of those dict contains hyperopt defined parameter spaces (e.g. ``hp.choice(..)``)
    which will be passed to the corresponding object : model, trainer or train method
    when performing hyperoptimization.
    :param gene_dataset: scVI gene dataset
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: Trainer class (e.g ``UnsupervisedTrainer``)
    :param model_specific_kwargs: dict of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: dict of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param use_batches: If False, pass n_batch=0 to model else pass gene_dataset.n_batches
    :param is_best_training: True if training the model with the best hyperparameters
    :return: best value of the early stopping metric, and best model if is_best_training
    """
    # hyperopt params
    space = defaultdict(dict, space)
    model_tunable_kwargs = space['model_tunable_kwargs']
    trainer_tunable_kwargs = space['trainer_tunable_kwargs']
    train_func_tunable_kwargs = space['train_func_tunable_kwargs']

    # add hardcoded parameters
    # disable scVI progbar
    trainer_specific_kwargs['show_progbar'] = False
    if is_best_training:
        trainer_specific_kwargs['train_size'] = 1.0
        # no monitoring, will crash otherwise
        trainer_specific_kwargs['frequency'] = None
        trainer_specific_kwargs['early_stopping_kwargs'] = None
    else:
        # evaluate at each epoch
        trainer_specific_kwargs['frequency'] = 1

    # merge params with fixed param precedence
    model_tunable_kwargs.update(model_specific_kwargs)
    trainer_tunable_kwargs.update(trainer_specific_kwargs)
    train_func_tunable_kwargs.update(train_func_specific_kwargs)

    if verbose and not is_best_training:
        print('Parameters being tested: ')
        print('model:')
        print(model_tunable_kwargs)
        print('trainer: ')
        print(trainer_tunable_kwargs)
        print('train method: ')
        print(train_func_tunable_kwargs)

    # define model
    model = model_class(
        n_input=gene_dataset.nb_genes,
        n_batch=gene_dataset.n_batches * use_batches,
        **model_tunable_kwargs,
    )

    # define trainer
    trainer = trainer_class(
        model,
        gene_dataset,
        **trainer_tunable_kwargs,
    )

    # train model
    trainer.train(**train_func_tunable_kwargs)

    # if training the best model, return model else return criterion
    if is_best_training:
        return trainer
    else:
        # select metric from early stopping kwargs or default to 'll_test_set'
        metric = None
        if 'early_stopping_kwargs' in trainer_specific_kwargs:
            early_stopping_kwargs = trainer_specific_kwargs['early_stopping_kwargs']
            if 'early_stopping_metric' in early_stopping_kwargs:
                metric = early_stopping_kwargs['early_stopping_metric']
                metric += '_' + trainer.early_stopping.on
        metric = metric if metric else 'll_test_set'

        return {
            'loss': trainer.history[metric][-1],
            'status': STATUS_OK,
            'history': trainer.history,
            'space': space,
        }
