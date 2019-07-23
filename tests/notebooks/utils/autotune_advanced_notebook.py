import copy
import datetime
import os
import pickle
from collections import defaultdict
from functools import partial
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from colour import Color
from hyperopt import STATUS_OK
from hyperopt import Trials
from scvi.inference import SemiSupervisedTrainer
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI
from scvi.models import VAE
from sklearn.model_selection import train_test_split


def custom_objective_hyperopt(
    space, is_best_training=False, dataset=None, n_epochs=None
):
    """Custom objective function for advanced autotune tutorial."""
    space = defaultdict(dict, space)
    model_tunable_kwargs = space["model_tunable_kwargs"]
    trainer_tunable_kwargs = space["trainer_tunable_kwargs"]
    train_func_tunable_kwargs = space["train_func_tunable_kwargs"]

    trainer_specific_kwargs = {}
    model_specific_kwargs = {}
    train_func_specific_kwargs = {}
    trainer_specific_kwargs["use_cuda"] = bool(torch.cuda.device_count())
    train_func_specific_kwargs["n_epochs"] = n_epochs

    # add hardcoded parameters
    # disable scVI progbar
    trainer_specific_kwargs["show_progbar"] = False
    trainer_specific_kwargs["frequency"] = 1

    # merge params with fixed param precedence
    model_tunable_kwargs.update(model_specific_kwargs)
    trainer_tunable_kwargs.update(trainer_specific_kwargs)
    train_func_tunable_kwargs.update(train_func_specific_kwargs)

    scanvi = SCANVI(
        dataset.nb_genes, dataset.n_batches, dataset.n_labels, **model_tunable_kwargs
    )
    trainer_scanvi = SemiSupervisedTrainer(scanvi, dataset, **trainer_tunable_kwargs)
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
        indices=np.squeeze(dataset.batch_indices == 1)
    )
    trainer_scanvi.unlabelled_set.to_monitor = ["reconstruction_error", "accuracy"]
    indices_labelled = np.squeeze(dataset.batch_indices == 0)

    if not is_best_training:
        # compute k-fold accuracy on a 20% validation set
        k = 5
        accuracies = np.zeros(k)
        indices_labelled = np.squeeze(dataset.batch_indices == 0)
        for i in range(k):
            indices_labelled_train, indices_labelled_val = train_test_split(
                indices_labelled.nonzero()[0], test_size=0.2
            )
            trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(
                indices=indices_labelled_train
            )
            trainer_scanvi.labelled_set.to_monitor = [
                "reconstruction_error",
                "accuracy",
            ]
            trainer_scanvi.validation_set = trainer_scanvi.create_posterior(
                indices=indices_labelled_val
            )
            trainer_scanvi.validation_set.to_monitor = ["accuracy"]
            trainer_scanvi.train(**train_func_tunable_kwargs)
            accuracies[i] = trainer_scanvi.history["accuracy_unlabelled_set"][-1]
        return {"loss": -accuracies.mean(), "space": space, "status": STATUS_OK}
    else:
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(
            indices=indices_labelled
        )
        trainer_scanvi.labelled_set.to_monitor = ["reconstruction_error", "accuracy"]
        trainer_scanvi.train(**train_func_tunable_kwargs)
        return trainer_scanvi


class Benchmarkable:
    """Helper class for Hyperparameter tuning blog post."""

    def __init__(
        self,
        global_path: str = None,
        exp_key: str = None,
        trainer_fname: str = None,
        model_fname: str = None,
        trials_fname: str = None,
        name: str = None,
        is_one_shot: bool = False,
        is_trials_only: bool = False,
    ):
        # indicator
        self.is_one_shot = is_one_shot
        self.is_trials_only = is_trials_only

        # construct paths if need be
        if global_path and exp_key and not trainer_fname and not is_trials_only:
            trainer_fname = os.path.join(global_path, "best_trainer_" + exp_key)
        if global_path and exp_key and not model_fname and not is_trials_only:
            model_fname = os.path.join(global_path, "best_model_" + exp_key)
        if global_path and exp_key and not trials_fname and not is_one_shot:
            trials_fname = os.path.join(global_path, "trials_" + exp_key)

        # load pickled attributes
        location = "cuda:0" if torch.cuda.is_available() else "cpu"
        with open(trainer_fname, "rb") as f:
            self.trainer: UnsupervisedTrainer = pickle.load(f)
        with open(model_fname, "rb") as f:
            self.model: VAE = torch.load(f, map_location=location)
        if not is_one_shot:
            with open(trials_fname, "rb") as f:
                self.trials: Trials = pickle.load(f)

        # runtime attributes
        if not self.is_one_shot:
            self.n_evals = len(self.trials.trials)
            # runtime
            self.total_train_time = datetime.timedelta(
                seconds=sum(
                    [trial["result"]["elapsed_time"] for trial in self.trials.trials]
                )
            )
            self.runtime = (
                self.trials.trials[-1]["book_time"]
                + datetime.timedelta(
                    seconds=self.trials.trials[-1]["result"]["elapsed_time"]
                )
                - self.trials.trials[0]["book_time"]
            )
            # elbo history
            self.history = self.trials.best_trial["result"]["history"]["elbo_test_set"]
            self.history_train = self.trials.best_trial["result"]["history"][
                "elbo_train_set"
            ]
            self.best_performance = self.trials.best_trial["result"]["loss"]
            self.n_epochs = self.trials.best_trial["result"]["space"][
                "train_func_tunable_kwargs"
            ].get("n_epochs", None)
            if not self.n_epochs:
                self.n_epochs = len(self.history)
        else:
            # elbo history is in trainer if one_shot
            self.history = self.trainer.history["elbo_test_set"]
            self.history_train = self.trainer.history["elbo_train_set"]
            self.best_performance = self.trainer.test_set.marginal_ll(n_mc_samples=100)
            self.n_epochs = len(self.history)

        self.name = name

        # imputed values
        self.original = None
        self.imputed = None

        # imputation scores
        self.median_score = None
        self.mean_score = None

    @property
    def imputation_errors(self):
        return np.abs(np.concatenate(self.original) - np.concatenate(self.imputed))

    def compute_imputed(
        self, n_epochs=None, n_samples=1, rate=0.1, corruption="uniform"
    ):
        # corrupt data with
        corrupted_trainer = copy.deepcopy(self.trainer)
        corrupted_trainer.train_test_validation(train_size=1.0)
        corrupted_trainer.corrupt_posteriors(rate=rate, corruption=corruption)
        corrupted_trainer.show_progbar = True
        # use n_epochs from early stopping
        if n_epochs is None:
            n_epochs = self.n_epochs
        corrupted_trainer.train(n_epochs=n_epochs)
        corrupted_trainer.uncorrupt_posteriors()

        (original_list, imputed_list) = corrupted_trainer.train_set.imputation_list(
            n_samples=n_samples
        )

        # Mean/Median of medians for each cell
        imputation_cells = []
        for original, imputed in zip(original_list, imputed_list):
            has_imputation = len(original) and len(imputed)
            imputation_cells += [
                np.median(np.abs(original - imputed)) if has_imputation else 0
            ]
        self.median_score = np.median(imputation_cells)
        self.mean_score = np.mean(imputation_cells)
        print(
            "\nMedian of Median: %.4f\nMean of Median for each cell: %.4f"
            % (self.median_score, self.mean_score)
        )

        # store imputed values
        self.imputed = imputed_list
        self.original = original_list

    def get_param_df(self):
        ddd = {}
        for i, trial in enumerate(self.trials):
            dd = {}
            dd["marginal_ll"] = trial["result"]["loss"]
            for item in trial["result"]["space"].values():
                for key, value in item.items():
                    dd[key] = value
            ddd[i] = dd
        df_space = pd.DataFrame(ddd)
        df_space = df_space.T
        n_params_dataset = np.vectorize(
            partial(n_params, self.trainer.gene_dataset.nb_genes)
        )
        df_space["n_params"] = n_params_dataset(
            df_space["n_layers"], df_space["n_hidden"], df_space["n_latent"]
        )
        df_space = df_space[
            [
                "marginal_ll",
                "n_layers",
                "n_hidden",
                "n_latent",
                "reconstruction_loss",
                "dropout_rate",
                "lr",
                "n_epochs",
                "n_params",
            ]
        ]
        df_space = df_space.sort_values(by="marginal_ll")
        df_space["run index"] = df_space.index
        df_space.index = np.arange(1, df_space.shape[0] + 1)
        return df_space


def n_params(n_genes, n_layers, n_hidden, n_latent):
    if n_layers == 0:
        res = 2 * n_genes * n_latent
    else:
        res = 2 * n_genes * n_hidden
        for i in range(n_layers - 1):
            res += 2 * n_hidden * n_hidden
        res += 2 * n_hidden * n_latent
    return res


class PlotBenchmarkables:
    """Helper class for Hyperparameter tuning blog post."""

    def __init__(
        self,
        tuned_benchmarkables: Dict[str, Benchmarkable],
        one_shot_benchmarkables: Dict[str, Benchmarkable],
    ):
        self.tuned_benchmarkables = tuned_benchmarkables
        self.one_shot_benchmarkables = one_shot_benchmarkables

    def get_runtime_dataframe(self):
        names = list(self.tuned_benchmarkables.keys())
        cols = [
            "Nb cells",
            "Nb genes",
            "Total GPU time",
            "Total wall time",
            "Number of trainings",
            "Avg time per training",
            "Avg epochs per training",
            "Number of GPUs",
            "Best epoch",
            "Max epoch",
        ]
        df_res = pd.DataFrame(index=names, columns=cols)
        for name in names:
            sub_df = df_res.loc[name]
            benchmarkable = self.tuned_benchmarkables[name]
            runtime_info = defaultdict(list)

            for trial in benchmarkable.trials:
                result = trial["result"]
                runtime_info["best_epoch"].append(result["best_epoch"])
                runtime_info["train_time"].append(result["elapsed_time"])

            def fill_sub_df(sub_df, benchmarkable):
                sub_df["Nb cells"] = benchmarkable.trainer.gene_dataset.X.shape[0]
                sub_df["Nb genes"] = benchmarkable.trainer.gene_dataset.nb_genes
                sub_df["Total GPU time"] = benchmarkable.total_train_time
                sub_df["Total wall time"] = benchmarkable.runtime
                sub_df["Best epoch"] = benchmarkable.n_epochs
                sub_df["Max epoch"] = (
                    50 if benchmarkable.name == "Brain Large tuned" else 1000
                )
                sub_df["Number of trainings"] = benchmarkable.n_evals
                sub_df["Avg time per training"] = np.mean(runtime_info["train_time"])
                sub_df["Avg epochs per training"] = np.mean(runtime_info["best_epoch"])
                sub_df["Number of GPUs"] = (
                    16 if benchmarkable.name == "Brain Large tuned" else 1
                )

            fill_sub_df(sub_df, benchmarkable)

        return df_res

    def get_results_dataframe(self):
        types = ["tuned", "default"]
        names = list(self.tuned_benchmarkables.keys())
        index = [(name, t) for name in names for t in types]
        index = pd.MultiIndex.from_tuples(index)

        cols = [
            ("Likelihood", t)
            for t in ["Held-out marginal ll", "ELBO train", "ELBO test"]
        ]
        cols += [("Imputation score", t) for t in ["median", "mean"]]
        cols = pd.MultiIndex.from_tuples(cols)

        df_res = pd.DataFrame(index=index, columns=cols)
        for name in names:
            tuned_benchmarkable = self.tuned_benchmarkables[name]
            one_shot_benchmarkable = self.one_shot_benchmarkables[name]

            def fill_sub_df(sub_df, benchmarkable):
                sub_df[("Imputation score", "median")] = benchmarkable.median_score
                sub_df[("Imputation score", "mean")] = benchmarkable.mean_score
                sub_df[
                    ("Likelihood", "Held-out marginal ll")
                ] = benchmarkable.best_performance
                sub_df[("Likelihood", "ELBO train")] = benchmarkable.history_train[
                    benchmarkable.n_epochs - 1
                ]
                sub_df[("Likelihood", "ELBO test")] = benchmarkable.history[
                    benchmarkable.n_epochs - 1
                ]

            fill_sub_df(df_res.loc[(name, "tuned")], tuned_benchmarkable)
            fill_sub_df(df_res.loc[(name, "default")], one_shot_benchmarkable)

        return df_res

    def plot_histories(
        self,
        ylims_dict=None,
        alpha=0.1,
        save_path="",
        filename="elbo_histories",
        **fig_kwargs
    ):
        fig, axes = plt.subplots(
            ncols=len(self.tuned_benchmarkables), nrows=1, **fig_kwargs
        )
        for i, name in enumerate(self.tuned_benchmarkables.keys()):
            tuned_benchmarkable = self.tuned_benchmarkables[name]
            one_shot_benchmarkable = self.one_shot_benchmarkables[name]
            self.plot_histories_single(
                tuned_benchmarkable, ax=axes[i], name=name, alpha=alpha
            )
            default_handles = axes[i].plot(
                one_shot_benchmarkable.trainer.history["elbo_test_set"],
                c="b",
                linewidth=1,
                label="default",
            )
            axes[i].legend(handles=default_handles)
            if ylims_dict:
                axes[i].set_ylim(ylims_dict[name])
        plt.savefig(os.path.join(save_path, filename))

    def plot_histories_single(
        self, benchmarkable, ax, name=None, ylim=None, xlim=None, alpha=0.1
    ):
        histories = []
        for trial in benchmarkable.trials.trials:
            histories.append(trial["result"]["history"]["elbo_test_set"])
        benchmarkable.history_df = pd.DataFrame(histories).T
        red = Color("red")
        colors = list(red.range_to(Color("green"), benchmarkable.n_evals))
        colors = [c.get_rgb() for c in colors]
        benchmarkable.history_df.plot(
            ax=ax,
            ylim=ylim,
            xlim=xlim,
            color=colors,
            alpha=alpha,
            legend=False,
            label=None,
        )
        name = name if name else benchmarkable.name
        ax.set_title(name + " : Held-out ELBO")
        ax.set_xlabel("epoch")
        ax.set_xlabel("ELBO")
