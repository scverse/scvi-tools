import numpy as np

import scvi
from scvi import REGISTRY_KEYS


class TestSemiSupervisedTrainingPlan(scvi.train.SemiSupervisedTrainingPlan):
    def __init__(self, *args, **kwargs):
        self.labels_to_indices = kwargs.pop("labels_to_indices")
        self.n_samples_per_label = kwargs.pop("n_samples_per_label")

        self.current_epoch_labeled_indices = None

        super().__init__(*args, **kwargs)

    def training_step(self, batch, batch_idx):
        _, labeled_tensors = batch
        labels, _ = (
            labeled_tensors[REGISTRY_KEYS.LABELS_KEY].cpu().numpy(),
            labeled_tensors[REGISTRY_KEYS.BATCH_KEY].cpu().numpy(),
        )
        if self.n_samples_per_label is None:
            return super().training_step(batch, batch_idx)

        _, counts = np.unique(labels, return_counts=True)
        assert np.all(counts == self.n_samples_per_label)

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
    labels_to_indices = {}
    for index, label in adata.obs[["indices", "labels"]].values:
        if label not in labels_to_indices:
            labels_to_indices[label] = []
        labels_to_indices[label].append(index)

    scvi.model.SCANVI._training_plan_cls = TestSemiSupervisedTrainingPlan
    plan_kwargs = {
        "n_samples_per_label": n_samples_per_label,
        "labels_to_indices": labels_to_indices,
    }
    scvi.model.SCANVI.setup_anndata(
        adata,
        batch_key="indices",  # hack to load indices
        labels_key="labels",
        unlabeled_category="label_0",
    )
    model = scvi.model.SCANVI(adata)
    model.train(
        max_epochs=2,
        batch_size=batch_size // 2,
        n_samples_per_label=n_samples_per_label,
        plan_kwargs=plan_kwargs,
    )
