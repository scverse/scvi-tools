import numpy as np
import pytest
import torch

import scvi
from scvi import REGISTRY_KEYS
from tests.dataset.utils import generic_setup_adata_manager


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


def test_anndataloader_distributed_sampler_init():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    with pytest.raises(ValueError):
        _ = scvi.dataloaders.AnnDataLoader(
            manager,
            sampler="a sampler",
            distributed_sampler=True,
        )


def test_contrastive_dataloader(
    mock_contrastive_adata_manager,
    mock_background_indices,
    mock_target_indices,
):
    batch_size = 32
    adata = mock_contrastive_adata_manager.adata

    expected_background_indices = mock_background_indices[:batch_size]
    expected_background_input = torch.Tensor(
        adata.layers["raw_counts"][expected_background_indices, :]
    )
    expected_background_labels = torch.LongTensor(
        adata.obs["_scvi_labels"].iloc[expected_background_indices]
    ).unsqueeze(1)
    expected_background_batch = torch.LongTensor(
        adata.obs["_scvi_batch"].iloc[expected_background_indices]
    ).unsqueeze(1)

    expected_target_indices = mock_target_indices[:batch_size]
    expected_target_input = torch.Tensor(
        adata.layers["raw_counts"][expected_target_indices, :]
    )
    expected_target_labels = torch.LongTensor(
        adata.obs["_scvi_labels"].iloc[expected_target_indices]
    ).unsqueeze(1)
    expected_target_batch = torch.LongTensor(
        adata.obs["_scvi_batch"].iloc[expected_target_indices]
    ).unsqueeze(1)

    dataloader = scvi.dataloaders.ContrastiveDataLoader(
        adata_manager=mock_contrastive_adata_manager,
        background_indices=mock_background_indices,
        target_indices=mock_target_indices,
        batch_size=batch_size,
    )
    batch = next(batch for batch in dataloader)

    assert isinstance(batch, dict)
    assert len(batch.keys()) == 2
    assert "background" in batch.keys()
    assert "target" in batch.keys()

    assert torch.equal(
        batch["background"][REGISTRY_KEYS.X_KEY], expected_background_input
    )
    assert torch.equal(
        batch["background"][REGISTRY_KEYS.LABELS_KEY], expected_background_labels
    )
    assert torch.equal(
        batch["background"][REGISTRY_KEYS.BATCH_KEY], expected_background_batch
    )

    assert torch.equal(batch["target"][REGISTRY_KEYS.X_KEY], expected_target_input)
    assert torch.equal(
        batch["target"][REGISTRY_KEYS.LABELS_KEY], expected_target_labels
    )
    assert torch.equal(batch["target"][REGISTRY_KEYS.BATCH_KEY], expected_target_batch)


def multiprocessing_worker(
    rank: int, world_size: int, manager: scvi.data.AnnDataManager, save_path: str
):
    # initializes the distributed backend that takes care of synchronizing processes
    torch.distributed.init_process_group(
        "gloo",  # backend that works on all systems
        init_method=f"file://{save_path}/dist_file",
        rank=rank,
        world_size=world_size,
    )

    _ = scvi.dataloaders.AnnDataLoader(manager, distributed_sampler=True)

    return


@pytest.mark.optional
def test_anndataloader_distributed_sampler(save_path: str, num_processes: int = 2):
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    torch.multiprocessing.spawn(
        multiprocessing_worker,
        args=(num_processes, manager, save_path),
        nprocs=num_processes,
        join=True,
    )
