# Optimize SCVI model with hyperparameter tuning

:::{note}
 to run scvi-tools with hyperparameters' tuning support, use: pip install scvi-tools[autotune]
:::

Hyperparameter tuning is the process of adjusting the parameters that control the training process of a machine learning model to find the best configuration for achieving optimal performance. These hyperparameters could include learning rate, batch size, the number of layers, and more. In PyTorch Lightning, when using Ray for hyperparameter tuning, you can leverage [Ray Tune](https://docs.ray.io/en/latest/tune/index.html), which is a scalable library for distributed hyperparameter optimization. To perform hyperparameter tuning in PyTorch Lightning with Ray, you first define a search space for the hyperparameters you want to tune (such as learning rate or batch size). Then, you set up a TuneReportCallback to track the performance of each training run and report the results back to Ray Tune. Ray will then automatically run multiple trials with different hyperparameter combinations and help you find the best-performing set.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/use_cases/autotune_scvi`
```

There are several common parameters that need to be entered when running hyperparameters tuning with ray:
- model_cls: Model class on which to tune hyperparameters.
- metrics: Either a single metric or a list of metrics to track during the experiment and select the best hyperparamters based on.
We typically use a metric that is one from the common metrics we calculate at each step, like "elbo_validation".
Starting scvi-tools v1.3, we are now able to search for optimal [scib-metrics](https://scib-metrics.readthedocs.io/en/stable/) of any choice (metrics that can qualify the quality of the batch correction and bio conservation results, or a mixture between them).
List of metrics names that can be used during scib-metrics autotune: [Total](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.benchmark.Benchmarker.html#scib_metrics.benchmark.Benchmarker) (everything, default), all [Batch correction](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.benchmark.BatchCorrection.html#scib_metrics.benchmark.BatchCorrection) metrics, all [Bio conservation](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.benchmark.BioConservation.html#scib_metrics.benchmark.BioConservation) metrics, [Silhouette label](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.silhouette_label.html#scib_metrics.silhouette_label), [Isolated labels](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.isolated_labels.html#scib_metrics.isolated_labels), [Leiden](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html#scib_metrics.nmi_ari_cluster_labels_leiden), [KMeans](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.nmi_ari_cluster_labels_kmeans.html#scib_metrics.nmi_ari_cluster_labels_kmeans), [cLISI](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.clisi_knn.html#scib_metrics.clisi_knn), [iLISI](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.ilisi_knn.html#scib_metrics.ilisi_knn), [KBET](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.kbet_per_label.html#scib_metrics.kbet), [BRAS](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.bras.html#scib_metrics.bras), [Graph connectivity](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.graph_connectivity.html#scib_metrics.graph_connectivity), [PCR comparison](https://scib-metrics.readthedocs.io/en/latest/generated/scib_metrics.pcr_comparison.html#scib_metrics.pcr_comparison)
- num_samples: Total number of hyperparameter configurations to sample from the search space.
- search_space: Dictionary of hyperparameter names and their respective search spaces. can contain any parameter from the model's parameter (parameters to pass to the model constructor) and any from the training parameters (parameters to pass to the model's ``train`` method)
- mode: one of ``'min'``, ``'max'``. Is a model preferable that minimizes or maximizes the metrics?
- data: the data that we are using and was setup with the model.

Example of running an autotune experiment for a SCVI model, with a search space of max_epochs and n_hidden:
```python
adata = synthetic_iid()
SCVI.setup_anndata(adata)

experiment = run_autotune(
    model_cls=SCVI,
    data=adata,
    metrics=["elbo_validation"],
    mode="min",
    search_space={
        "model_params": {
            "n_hidden": tune.choice([1, 2]),
        },
        "train_params": {
            "max_epochs": tune.choice([1, 20]),
        },
    },
    num_samples=2,
    seed=0,
    scheduler="asha",
    searcher="hyperopt",
)
```
For a more end to end example, see our scvi autotune tutorial {doc}`/tutorials/notebooks/use_cases/autotune_scvi`
