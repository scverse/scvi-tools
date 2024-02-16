from __future__ import annotations

from typing import Literal

import lightning.pytorch as pl
from anndata import AnnData
from lightning.pytorch.callbacks import Callback
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation


class ScibCallback(Callback):
    def __init__(
        self,
        bio_conservation_metrics: BioConservation | None = None,
        batch_correction_metrics: BatchCorrection | None = None,
        stage: Literal["training", "validation", "both"] = "both",
    ):
        super().__init__()

        self.bio_conservation_metrics = bio_conservation_metrics
        self.batch_correction_metrics = batch_correction_metrics
        self.stage = stage

    def compute_metrics(
        self,
        trainer: pl.Trainer,
        pl_module: pl.LightningModule,
        stage: Literal["training", "validation"],
    ):
        if stage == "training" and self.stage not in ["training", "both"]:
            return
        elif stage == "validation" and self.stage not in ["validation", "both"]:
            return

        if not hasattr(pl_module, f"_{stage}_epoch_outputs"):
            raise ValueError(f"The training plan must have a `_{stage}_epoch_outputs` attribute.")

        outputs = getattr(pl_module, f"_{stage}_epoch_outputs")
        x = outputs["x"].numpy()
        z = outputs["z"].numpy()
        batch = outputs["batch"].numpy()
        labels = outputs["labels"].numpy()

        adata = AnnData(X=x, obs={"batch": batch, "labels": labels}, obsm={"z": z})
        benchmarker = Benchmarker(
            adata,
            batch_key="batch",
            label_key="labels",
            embedding_obsm_keys=["z"],
            bio_conservation_metrics=self.bio_conservation_metrics,
            batch_correction_metrics=self.batch_correction_metrics,
        )
        benchmarker.prepare()
        benchmarker.benchmark()
        results = benchmarker.get_results(min_max_scale=False).to_dict()
        metrics = {f"training {metric}": results[metric]["z"] for metric in results}
        pl_module.logger.log_metrics(metrics, trainer.global_step)

        delattr(pl_module, f"_{stage}_epoch_outputs")

    def on_train_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self.compute_metrics(trainer, pl_module, "training")

    def on_validation_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self.compute_metrics(trainer, pl_module, "validation")
