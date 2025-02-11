import pytest
from ray import tune
from ray.tune import ResultGrid

from scvi import settings
from scvi.autotune import AutotuneExperiment, run_autotune
from scvi.data import synthetic_iid
from scvi.dataloaders import DataSplitter
from scvi.model import SCVI

# n_batches=3
# adata = synthetic_iid()
# SCVI.setup_anndata(adata, batch_key="batch")
# model = SCVI(adata,n_latent=5)
# model.train(max_epochs=1)
# SCVI_LATENT_KEY = "X_scVI"
# adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
# bm = Benchmarker(
#     adata,
#     batch_key="batch",
#     label_key="labels",
#     embedding_obsm_keys=[SCVI_LATENT_KEY],
#     #n_jobs=-1,
# )
# bm.prepare()
# bm.benchmark()
# results = bm.get_results(min_max_scale=False).to_dict()
# metrics = {f"training {metric}": results[metric]["z"] for metric in results}


def test_run_autotune_scvi_basic(save_path: str):
    settings.logging_dir = save_path
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)

    experiment = run_autotune(
        SCVI,
        adata,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


def test_run_autotune_scvi_no_anndata(save_path: str, n_batches: int = 3):
    settings.logging_dir = save_path
    adata = synthetic_iid(n_batches=n_batches)
    SCVI.setup_anndata(adata, batch_key="batch")
    manager = SCVI._get_most_recent_anndata_manager(adata)

    datamodule = DataSplitter(manager)
    datamodule.n_vars = adata.n_vars
    datamodule.n_batch = n_batches

    experiment = run_autotune(
        SCVI,
        data=datamodule,
        metrics=["elbo_validation"],
        mode="min",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


@pytest.mark.parametrize("metric", ["Total", "Bio conservation", "iLISI"])
def test_run_autotune_scvi_with_scib(metric: str, save_path: str = "."):
    settings.logging_dir = save_path
    adata = synthetic_iid(batch_size=10, n_genes=10)
    SCVI.setup_anndata(adata)

    # Mapping of metric fn names to clean DataFrame column names
    # metric_name_cleaner = {
    #     "silhouette_label": "Silhouette label",
    #     "silhouette_batch": "Silhouette batch",
    #     "isolated_labels": "Isolated labels",
    #     "nmi_ari_cluster_labels_leiden_nmi": "Leiden NMI", - not run
    #     "nmi_ari_cluster_labels_leiden_ari": "Leiden ARI", - not run
    #     "nmi_ari_cluster_labels_kmeans_nmi": "KMeans NMI",
    #     "nmi_ari_cluster_labels_kmeans_ari": "KMeans ARI",
    #     "clisi_knn": "cLISI",
    #     "ilisi_knn": "iLISI",
    #     "kbet_per_label": "KBET",
    #     "graph_connectivity": "Graph connectivity",
    #     "pcr_comparison": "PCR comparison",
    #     "BatchCorrection": "Batch correction",
    #     "BioConservation": "Bio conservation",
    #     "SCIB_Total": "Total",
    # }

    experiment = run_autotune(
        SCVI,
        adata,
        metrics=[metric],
        mode="max",
        search_space={
            "model_params": {
                "n_hidden": tune.choice([1, 2]),
            },
            "train_params": {
                "max_epochs": 1,
            },
        },
        num_samples=2,
        seed=0,
        scheduler="asha",
        searcher="hyperopt",
        scib_stage="validation",
        scib_subsample_rows=50,
    )
    assert isinstance(experiment, AutotuneExperiment)
    assert hasattr(experiment, "result_grid")
    assert isinstance(experiment.result_grid, ResultGrid)


# def test_early_stopping():
#     # we use this temporarily to debug the scib-metrics callback
#     n_epochs = 100
#
#     adata = synthetic_iid()
#     SCVI.setup_anndata(
#         adata,
#         batch_key="batch",
#         labels_key="labels",
#     )
#     model = SCVI(adata)
#     model.train(n_epochs, early_stopping=True, plan_kwargs={"lr": 0})
#     assert len(model.history["elbo_train"]) < n_epochs
