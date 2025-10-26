import os
import pickle

import pandas as pd
import torch

from scvi.utils import dependencies


def _compute_kl_weight(
    epoch: int,
    step: int,
    n_epochs_kl_warmup: int | None,
    n_steps_kl_warmup: int | None,
    max_kl_weight: float = 1.0,
    min_kl_weight: float = 0.0,
) -> float | torch.Tensor:
    """Computes the kl weight for the current step or epoch.

    If both `n_epochs_kl_warmup` and `n_steps_kl_warmup` are None `max_kl_weight` is returned.

    Parameters
    ----------
    epoch
        Current epoch.
    step
        Current step.
    n_epochs_kl_warmup
        Number of training epochs to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`
    n_steps_kl_warmup
        Number of training steps (minibatches) to scale weight on KL divergences from
        `min_kl_weight` to `max_kl_weight`
    max_kl_weight
        Maximum scaling factor on KL divergence during training.
    min_kl_weight
        Minimum scaling factor on KL divergence during training.
    """
    if min_kl_weight > max_kl_weight:
        raise ValueError(
            f"min_kl_weight={min_kl_weight} is larger than max_kl_weight={max_kl_weight}."
        )

    slope = max_kl_weight - min_kl_weight
    if n_epochs_kl_warmup:
        if epoch < n_epochs_kl_warmup:
            return slope * (epoch / n_epochs_kl_warmup) + min_kl_weight
    elif n_steps_kl_warmup:
        if step < n_steps_kl_warmup:
            return slope * (step / n_steps_kl_warmup) + min_kl_weight
    return max_kl_weight


def _safe_load_logger_history(trainer):
    hist = getattr(trainer.logger, "history", None)
    if hist:
        return {k: v.copy() for k, v in hist.items()}  # deep copy from memory
    history_path = getattr(trainer.logger, "history_path", None)  #  file (written by rank-0)
    if history_path and os.path.exists(history_path):
        with open(history_path, "rb") as f:
            return pickle.load(f)
    return None


@dependencies("mlflow")
def _mlflow_logger(
    model=None, trainer=None, training_plan=None, data_splitter=None, run_id: str = None
):
    import mlflow

    # Log training/architecture parameters
    if trainer is not None:
        mlflow.log_params({"max_epochs": trainer.max_epochs}, run_id=run_id)
    if training_plan is not None:
        mlflow.log_params({"lr": training_plan.lr}, run_id=run_id)
    if data_splitter is not None:
        mlflow.log_params(data_splitter.data_loader_kwargs, run_id=run_id)
        mlflow.log_params(data_splitter.data_loader_kwargs, run_id=run_id)
        mlflow.log_params(data_splitter.data_loader_kwargs, run_id=run_id)

    if model is not None:
        mlflow.log_params(model.init_params_["non_kwargs"], run_id=run_id)
        mlflow.log_params(model.init_params_["kwargs"], run_id=run_id)
        mlflow.log_params(model.registry["setup_args"], run_id=run_id)

        #     {
        #         # "train_size": train_size,
        #         # "check_val_every_n_epoch": check_val_every_n_epoch,
        #     }
        # )

        # Log experiment metrics (from model.history)
        for key, values in model.history.items():
            if isinstance(values, (list, tuple)):
                for i, v in enumerate(values):
                    mlflow.log_metric(key, v, step=i, run_id=run_id)
            # If it's a pandas Series
            elif isinstance(values, pd.Series):
                for i, v in enumerate(values):
                    if pd.notna(v):
                        mlflow.log_metric(key, float(v), step=i, run_id=run_id)
            # Pandas DataFrame
            elif isinstance(values, pd.DataFrame):
                # Log all numeric columns as metrics
                for col in values.columns:
                    for i, v in enumerate(values[col]):
                        if pd.notna(v):
                            mlflow.log_metric(col, float(v), step=i, run_id=run_id)
            # If it's a scalar
            else:
                try:
                    mlflow.log_metric(key, float(values), run_id=run_id)
                except (TypeError, ValueError):
                    # Non-numeric scalar, store as param instead
                    mlflow.log_param(key, str(values), run_id=run_id)
    return
