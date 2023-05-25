import numpy as np

import scvi
from scvi import REGISTRY_KEYS


class TestSemiSupervisedTrainingPlan(scvi.train.SemiSupervisedTrainingPlan):
    def __init__(self, *args, **kwargs):
        self.n_samples_per_label = kwargs.pop("n_samples_per_label")
        self.epoch_to_labeled_indices = {}

        super().__init__(*args, **kwargs)

    def training_step(self, batch, batch_idx):
        _, labeled_tensors = batch
        labels, indices = (
            labeled_tensors[REGISTRY_KEYS.LABELS_KEY].cpu().numpy(),
            labeled_tensors[REGISTRY_KEYS.BATCH_KEY].cpu().numpy(),
        )
        if self.n_samples_per_label is None:
            return super().training_step(batch, batch_idx)

        _, counts = np.unique(labels, return_counts=True)
        assert np.all(counts == self.n_samples_per_label)

        # new epoch
        if self.current_epoch not in self.epoch_to_labeled_indices:
            self.epoch_to_labeled_indices[self.current_epoch] = []

            prev_epoch = self.current_epoch - 1
            prev_prev_epoch = self.current_epoch - 2

            # check previous two epochs have different labeled indices
            if (
                prev_epoch in self.epoch_to_labeled_indices
                and prev_prev_epoch in self.epoch_to_labeled_indices
            ):
                prev_indices = self.epoch_to_labeled_indices[prev_epoch]
                prev_prev_indices = self.epoch_to_labeled_indices[prev_prev_epoch]

                assert len(np.setdiff1d(prev_indices, prev_prev_indices)) > 0

        self.epoch_to_labeled_indices[self.current_epoch] += indices.squeeze().tolist()

        return super().training_step(batch, batch_idx)


def test_semisuperviseddataloader_subsampling(
    batch_size: int = 128,
    n_batches: int = 2,
    n_labels: int = 3,
    n_samples_per_label: int = 10,
):
    adata = scvi.data.synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_labels=n_labels
    )
    adata.obs["indices"] = np.arange(adata.n_obs)

    original_training_plan_cls = scvi.model.SCANVI._training_plan_cls
    scvi.model.SCANVI._training_plan_cls = TestSemiSupervisedTrainingPlan
    plan_kwargs = {
        "n_samples_per_label": n_samples_per_label,
    }
    scvi.model.SCANVI.setup_anndata(
        adata,
        batch_key="indices",  # hack to load indices
        labels_key="labels",
        unlabeled_category="label_0",
    )
    model = scvi.model.SCANVI(adata)
    model.train(
        max_epochs=10,
        batch_size=batch_size // 2,
        n_samples_per_label=n_samples_per_label,
        plan_kwargs=plan_kwargs,
    )

    scvi.model.SCANVI._training_plan_cls = original_training_plan_cls
