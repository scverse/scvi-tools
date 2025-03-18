import os

import numpy as np
import pytest
import torch
from tests.data.utils import generic_setup_adata_manager

import scvi
from scvi import REGISTRY_KEYS
from scvi.model import SCANVI


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
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches, n_labels=n_labels)
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


def multiprocessing_worker(
    rank: int,
    world_size: int,
    manager: scvi.data.AnnDataManager,
    save_path: str,
    datasplitter_kwargs,
):
    # initializes the distributed backend that takes care of synchronizing processes
    torch.distributed.init_process_group(
        "nccl",
        init_method=f"file://{save_path}/dist_file",
        rank=rank,
        world_size=world_size,
        store=None,
    )

    _ = scvi.dataloaders.AnnDataLoader(manager, **datasplitter_kwargs)

    return


@pytest.mark.multigpu
@pytest.mark.parametrize("num_processes", [1, 2])
def test_anndataloader_distributed_sampler(num_processes: int, save_path: str):
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    file_path = save_path + "/dist_file"
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)

    torch.multiprocessing.spawn(
        multiprocessing_worker,
        args=(num_processes, manager, save_path, {}),
        nprocs=num_processes,
        join=True,
    )


@pytest.mark.multigpu
@pytest.mark.parametrize("num_processes", [1, 2])
def test_scanvi_with_distributed_sampler(num_processes: int, save_path: str):
    adata = scvi.data.synthetic_iid()
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    file_path = save_path + "/dist_file"
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)
    datasplitter_kwargs = {}
    # Multi-GPU settings
    datasplitter_kwargs["distributed_sampler"] = True
    datasplitter_kwargs["drop_last"] = False
    if num_processes == 1:
        datasplitter_kwargs["distributed_sampler"] = False
    model = SCANVI(adata, n_latent=10)

    # initializes the distributed backend that takes care of synchronizing processes
    torch.distributed.init_process_group(
        "nccl",  # backend that works on all systems
        init_method=f"file://{save_path}/dist_file",
        rank=0,
        world_size=num_processes,
        store=None,
    )

    model.train(1, datasplitter_kwargs=datasplitter_kwargs)

    torch.distributed.destroy_process_group()
