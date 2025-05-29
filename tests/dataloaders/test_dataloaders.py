import os

import numpy as np
import pytest
import torch
from anndata.experimental import AnnCollection
from scipy.sparse import csr_matrix
from tests.data.utils import generic_setup_adata_manager

import scvi
from scvi import REGISTRY_KEYS
from scvi.dataloaders import CollectionAdapter
from scvi.model import SCANVI, SCVI


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

    original_training_plan_cls = SCANVI._training_plan_cls
    SCANVI._training_plan_cls = TestSemiSupervisedTrainingPlan
    plan_kwargs = {
        "n_samples_per_label": n_samples_per_label,
    }
    SCANVI.setup_anndata(
        adata,
        batch_key="indices",  # hack to load indices
        labels_key="labels",
        unlabeled_category="label_0",
    )
    model = SCANVI(adata)
    model.train(
        max_epochs=10,
        batch_size=batch_size // 2,
        n_samples_per_label=n_samples_per_label,
        plan_kwargs=plan_kwargs,
    )

    SCANVI._training_plan_cls = original_training_plan_cls


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


def test_anncollection(save_path: str):
    adata1 = scvi.data.synthetic_iid()
    adata1.X = csr_matrix(adata1.X)

    adata2 = scvi.data.synthetic_iid()
    adata2.X = csr_matrix(adata2.X)

    # create a fake counts layer to test training
    adata1.layers["test"] = adata1.X
    adata2.layers["test"] = adata2.X

    # take annotations from the `pbmc` dataset and leave
    # annotations in `covid` as an Unknown
    adata2.obs["labels"] = "Unknown"

    # X is all raw counts
    assert np.all(np.mod(adata1.X[:10].toarray(), 1) == 0)
    assert np.all(np.mod(adata2.X[:10].toarray(), 1) == 0)

    # create an AnnCollection on a subset of the data
    adata = AnnCollection(
        [adata1, adata2],
        join_vars="inner",
        join_obs="inner",
        label="dataset",
    )
    print(adata)

    collection_adapter = CollectionAdapter(adata)

    SCVI.setup_anndata(
        collection_adapter,
        layer="test",
        batch_key="dataset",
        labels_key="labels",
    )

    model = SCVI(collection_adapter, n_latent=10)

    # we're only training for a few epochs to show it works
    model.train(max_epochs=1, check_val_every_n_epoch=1, train_size=0.9, early_stopping=True)

    # Generate cell representations # TODO: Not worknig yet
    # latent = model.get_latent_representation()
    # collection_adapter.obsm["X_scVI"] = latent

    # Save model (note we DONT save anndatas - as they are already on disk))
    dir_path = save_path + "/model_scvi_anncollection"
    model.save(dir_path, save_anndata=False, overwrite=True)

    # Load model again
    loaded_model = SCVI.load(dir_path, adata=collection_adapter)

    # create and prepare query data (TODO: need to support not same num of gene  for anncollection)
    adata3 = scvi.data.synthetic_iid()
    adata3.X = csr_matrix(adata3.X)
    adata4 = scvi.data.synthetic_iid()
    adata4.X = csr_matrix(adata4.X)

    # create a fake counts layer to test training
    adata3.layers["test"] = adata3.X
    adata4.layers["test"] = adata4.X

    # take annotations from the `pbmc` dataset and leave
    # annotations in `covid` as an Unknown
    adata4.obs["labels"] = "Unknown"

    # create an AnnCollection on a subset of the data
    adata_query = AnnCollection(
        [adata3, adata4],
        join_vars="inner",
        join_obs="inner",
        label="dataset",
    )
    print(adata_query)

    collection_adapter_query = CollectionAdapter(adata_query)

    SCVI.prepare_query_anndata(collection_adapter_query, dir_path)
    scvi_query = SCVI.load_query_data(collection_adapter_query, dir_path)
    scvi_query.train(max_epochs=1, plan_kwargs={"weight_decay": 0.0})

    # create a scanvi model from scvi model
    scanvi_model = SCANVI.from_scvi_model(
        loaded_model, "Unknown", labels_key="labels", adata=collection_adapter_query
    )
    scanvi_model.train(1)

    # create from the direct method
    SCANVI.setup_anndata(
        collection_adapter_query,
        layer="test",
        batch_key="dataset",
        labels_key="labels",
        unlabeled_category="Unknown",
    )

    model = scvi.model.SCANVI(collection_adapter_query, n_latent=10)

    # we're only training for a few epochs to show it works
    model.train(max_epochs=1)

    # TODO: Not worknig yet
    # latent = model.get_latent_representation()
    # collection_adapter_query.obsm["X_scanVI"] = latent

    predictions, attributions = model.predict(collection_adapter, ig_interpretability=True)
    ig_top_features_samples = attributions.head(5)
    print(ig_top_features_samples)
