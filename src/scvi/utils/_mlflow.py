import os
import sys
import warnings
from typing import Any, Union

import pandas as pd

from scvi import settings
from scvi.utils import dependencies


@dependencies("mlflow")
def mlflow_log_artifact(
    local_path: str,
    artifact_path: str | None = None,
    run_id: str | None = None,
    max_size_mb: float = 5.0,
) -> None:
    import mlflow

    if not os.path.isfile(local_path):
        raise FileNotFoundError(f"File not found: {local_path}")

    file_size_mb = os.path.getsize(local_path) / (1024 * 1024)
    if file_size_mb <= max_size_mb:
        mlflow.log_artifact(local_path, artifact_path=artifact_path, run_id=run_id)
    else:
        warnings.warn(
            f"File too large to log to MLFlow: {local_path}",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
    return


@dependencies("mlflow")
def mlflow_log_table(
    data: Union[dict[str, Any], "pd.DataFrame"],
    artifact_file: str | None = None,
    run_id: str | None = None,
    max_size_mb: float = 1.0,
) -> None:
    import mlflow

    if isinstance(data, pd.DataFrame):
        data["table_index"] = data.index
    file_size_mb = sys.getsizeof(data) / (1024 * 1024)
    if file_size_mb <= max_size_mb:
        mlflow.log_table(data, artifact_file=artifact_file, run_id=run_id)
    else:
        warnings.warn(
            "Object too large to log to MLFlow",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
    if isinstance(data, pd.DataFrame):
        del data["table_index"]
    return


@dependencies("mlflow")
def mlflow_log_text(
    text: str,
    artifact_file: str | None = None,
    run_id: str | None = None,
    max_size_mb: float = 1.0,
) -> None:
    import mlflow

    file_size_mb = sys.getsizeof(text) / (1024 * 1024)
    if file_size_mb <= max_size_mb:
        mlflow.log_text(text, artifact_file=artifact_file, run_id=run_id)
    else:
        warnings.warn(
            "Object too large to log to MLFlow",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
    return


@dependencies("mlflow")
def mlflow_logger(
    model=None, trainer=None, training_plan=None, data_splitter=None, run_id: str = None
):
    import mlflow

    # Log training/architecture parameters
    if trainer is not None:
        mlflow.log_params({"max_epochs": trainer.max_epochs}, run_id=run_id)
        mlflow.log_params(
            {"check_val_every_n_epoch": trainer.check_val_every_n_epoch}, run_id=run_id
        )
        mlflow.log_params({"log_every_n_steps": trainer.log_every_n_steps}, run_id=run_id)
        mlflow.log_params({"world_size": trainer.world_size}, run_id=run_id)
    if training_plan is not None:
        mlflow.log_params({"lr": training_plan.lr}, run_id=run_id)
        mlflow.log_params({"eps": training_plan.eps}, run_id=run_id)
        mlflow.log_params({"lr_factor": training_plan.lr_factor}, run_id=run_id)
        mlflow.log_params({"lr_patience": training_plan.lr_patience}, run_id=run_id)
        mlflow.log_params({"n_epochs_kl_warmup": training_plan.n_epochs_kl_warmup}, run_id=run_id)
        mlflow.log_params({"n_steps_kl_warmup": training_plan.n_steps_kl_warmup}, run_id=run_id)
        mlflow.log_params({"weight_decay": training_plan.weight_decay}, run_id=run_id)
        mlflow.log_params({"max_kl_weight": training_plan.max_kl_weight}, run_id=run_id)
        mlflow.log_params({"min_kl_weight": training_plan.min_kl_weight}, run_id=run_id)
    if data_splitter is not None:
        mlflow.log_params(data_splitter.data_loader_kwargs, run_id=run_id)
        mlflow.log_params({"train_size": data_splitter.train_size}, run_id=run_id)
        mlflow.log_params({"n_train": data_splitter.n_train}, run_id=run_id)
        mlflow.log_params({"n_val": data_splitter.n_val}, run_id=run_id)
    if model is not None:
        mlflow.log_params(model.init_params_.get("non_kwargs", {}), run_id=run_id)
        mlflow.log_params(model.init_params_.get("kwargs", {}), run_id=run_id)
        mlflow.log_params(model.registry.get("setup_args", {}), run_id=run_id)
        for field in model.registry["field_registries"].keys():
            mlflow.log_params(
                model.registry["field_registries"][field].get("summary_stats", {}), run_id=run_id
            )

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
