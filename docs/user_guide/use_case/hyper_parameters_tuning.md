# Optimize SCVI model with hyperparameter tuning

:::{note}
In order to run scvi-tools with hyperparameters tuning support, use: pip install scvi-tools[autotune]
:::

Hyperparameter tuning is the process of adjusting the parameters that control the training process of a machine learning model to find the best configuration for achieving optimal performance. These hyperparameters could include learning rate, batch size, the number of layers, and more. In SCVI-Tools we practically have mechanisms to optimize models:

## RAY

In PyTorch Lightning, when using Ray for hyperparameter tuning, you can leverage [Ray Tune](https://docs.ray.io/en/latest/tune/index.html), which is a scalable library for distributed hyperparameter optimization. To perform hyperparameter tuning in PyTorch Lightning with Ray, you first define a search space for the hyperparameters you want to tune (such as learning rate or batch size). Then, you set up a TuneReportCallback to track the performance of each training run and report the results back to Ray Tune. Ray will then automatically run multiple trials with different hyperparameter combinations and help you find the best-performing set.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/use_cases/autotune_scvi`
```

There are several common parameters that need to be entered when running hyper parameters tuning with ray:
- model_cls: Model class on which to tune hyperparameters.
- metrics: Either a single metric or a list of metrics to track during the experiment and select the best hyperparamters based on.
We typically use a metric that is one from the common metrics we calculate at each step, like "elbo_validation".
Starting scvi-tools v1.3, we are now able to search for optimal [scib-metrics](https://scib-metrics.readthedocs.io/en/stable/) of any choice (metrics that can qualify the quality of the batch correction and bio conservation results, or a mixture between them). This currently works only for SCVI and SCANVI models, while other models will be implemented based on user adoption.
- num_samples: Total number of hyperparameter configurations to sample from the search space.
- search_space: Dictionary of hyperparameter names and their respective search spaces. can contain any parameter from the model's parameter (parameters to pass to the model constructor) and any from the trainining parameters (parameters to pass to the model's ``train`` method)
- mode: one of ``'min'``, ``'max'``. Is a model preferable that minimizes or maximizes the metrics?
- data: the data that we are using and was setup with the model.

Example of running an autotune experiement for a SCVI model, with a search space of max_epochs and n_hidden:
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
For a more end to end example, see our scvi autotune tutorial (link above).

## MLFLOW

Starting scvi-tools v1.4.1, users can log their experiments in [mlflow](https://mlflow.org/), a popular MLOPS UI platform.
In order to use it, user should install mlflow sw as well as scvi-tools dependency:

```bash
pip install mlflow

pip install -U scvi-tools[mlflow]
```

Run mlflow on its server:

```bash
mlflow ui --host 0.0.0.0 --port 5000
```

Then supply the URL of its UI to scvi-tools settings at the beginning of code and you good to go.
The experiment will be logged with all possible metrics (model, training, systems and so on) and be compared to previous experiments.

For more information, See this [link](https://yosef-lab.notion.site/MLFLOW-SCVI-Tools-29957932739c809c94c1d866e64ca1cd)
