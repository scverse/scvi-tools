from __future__ import annotations

import lightning.pytorch as pl
from lightning.pytorch.callbacks import Callback
from scib_metrics.benchmark import BatchCorrection, BioConservation

from scvi import METRIC_KEYS, REGISTRY_KEYS
from scvi.utils import dependencies


class ScibCallback(Callback):
    @dependencies("scib_metrics")
    def __init__(
        self,
        batch_correction: BatchCorrection | None = None,
        bio_conservation: BioConservation | None = None,
    ):
        super().__init__()
        self.batch_correction = batch_correction
        self.bio_conservation = bio_conservation

    def on_validation_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        from anndata import AnnData
        from scib_metrics.benchmark import Benchmarker

        attr_key = "_validation_outputs"
        if not hasattr(pl_module, attr_key):
            raise ValueError(
                f"`ScibCallback` requires the training plan to set a `{attr_key}` "
                "attribute at the end of the validation step."
            )

        outputs = getattr(pl_module, attr_key)
        required_keys = [
            REGISTRY_KEYS.X_KEY,
            REGISTRY_KEYS.BATCH_KEY,
            REGISTRY_KEYS.LABELS_KEY,
            METRIC_KEYS.LATENT_KEY,
        ]
        if not all(key in outputs for key in required_keys):
            raise ValueError(
                "`ScibCallback` requires the validation outputs to contain the following keys: "
                f"{required_keys}."
            )

        adata = AnnData(
            X=outputs[REGISTRY_KEYS.X_KEY].cpu().numpy(),
            obs={
                REGISTRY_KEYS.BATCH_KEY: outputs[REGISTRY_KEYS.BATCH_KEY].cpu().numpy(),
                REGISTRY_KEYS.LABELS_KEY: outputs[REGISTRY_KEYS.LABELS_KEY].cpu().numpy(),
            },
            obsm={METRIC_KEYS.LATENT_KEY: outputs[METRIC_KEYS.LATENT_KEY].cpu().numpy()},
        )
        benchmarker = Benchmarker(
            adata,
            batch_key=REGISTRY_KEYS.BATCH_KEY,
            label_key=REGISTRY_KEYS.LABELS_KEY,
            embedding_obsm_keys=[METRIC_KEYS.LATENT_KEY],
            bio_conservation_metrics=self.bio_conservation,
            batch_correction_metrics=self.batch_correction,
        )
        benchmarker.prepare()
        benchmarker.benchmark()
        results = benchmarker.get_results(min_max_scale=False).to_dict()
        metrics = {
            metric.lower().replace(" ", "_"): results[metric][METRIC_KEYS.LATENT_KEY]
            for metric in results
        }
        pl_module.logger.log_metrics(metrics, step=trainer.global_step)

        delattr(pl_module, attr_key)
