# Train SCVI model with multi-GPU support

:::{note}
In order to run scvi-tools with mulyi-GPU support, use: pip install scvi-tools[cuda]
:::

SCVI-Tools v1.3.0 now supports training on a **multi GPU system**, which can significantly speed up training and allow you to handle larger datasets. It is supported only on Nvidia GPUs and DDP with CUDA backend.

CUDA (Compute Unified Device Architecture) is used to enable high-performance parallel computing on NVIDIA GPUs. DDP (Distributed Data Parallel) is a popular backend for training models on multiple GPUs simultaneously.

Using **CUDA** with **DDP** for multi-GPU training can significantly boost the performance and scalability of your model training, especially with large datasets and models. However, it also comes with its own set of challenges, such as the need for proper configuration, higher hardware costs, and potential software compatibility issues.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/use_cases/multiGPU`
```

## Some of the benefits of using MultiGPU training

### 1. **Faster Training**
   - **Parallel Processing:** Using multiple GPUs allows you to parallelize the training process, meaning different GPUs can work on different batches of data simultaneously, leading to faster training times.
   - **Scalability:** As you add more GPUs to the system, the training can scale up, further reducing the time it takes to train large models.
   - **Larger Models and Datasets:** With multiple GPUs, you can train models that are too large to fit into a single GPU’s memory. You can also work with larger datasets more efficiently.

See the following figure which shows the training time (y-axis) of a PBMC data, SCVI model with and without the use of multiGPU vs. different sizes of cells in different adata (x-axis), going from smallest to largest:
We can see the advantage with larger data, while for the small data, there's no much advantage due to multiGPU overhead needed.

:::{figure} /\_static/multigpu.png
:align: center
:alt: MultiGPU train time
:class: img-fluid
:::

### 2. **Efficiency with DDP**
   - **Data Parallelism:** DDP splits your data into smaller chunks, distributes them across multiple GPUs, and then aggregates the results. This allows for better use of available hardware resources.
   - **Synchronization:** DDP efficiently synchronizes gradients and updates between GPUs after each batch, ensuring that each GPU works in sync with the others without much manual intervention.
   - **Automatic Handling of Distributed Work:** DDP abstracts much of the complexity of managing multiple GPUs, so you don’t have to manually manage how data is split or how results are aggregated.

### 3. **Memory Utilization**
   - **Larger Memory Pool:** When using multiple GPUs, each GPU can hold a part of the model and data, effectively creating a larger memory pool. This allows for larger batch sizes or more complex models.
   - **Shared Memory Deduplication:** In subprocess-based DDP, each rank loads its own copy of `adata` into memory, which can be costly for large datasets. scvi-tools automatically uses POSIX shared memory to deduplicate `adata.X` across all ranks on the same node, so only one physical copy is kept in memory. For example, with 100K obs x 1000 genes (~400 MB), 8 GPUs go from 3.2 GB total to just 400 MB. This is auto-enabled when DDP is detected and works with both dense numpy and scipy sparse matrices.

## Using MultiGPU training in SCVI-Tools

1. Currently, multiGPU training is supported for the following models:
   - {class}`~scvi.model.SCVI`
   - {class}`~scvi.model.SCANVI`,
   - {class}`~scvi.model.CondSCVI`
   - {class}`~scvi.model.LinearSCVI`
   - {class}`~scvi.model.TOTALVI`,
   - {class}`~scvi.model.MULTIVI`
   - {class}`~scvi.model.PEAKVI`



2. How to use it: during the model train command, we need to use the strategy parameter, and its value depends on whether we are running from command line/script or from an interactive session like jupyter notebook/colab:
- From non-interactive session:
```python
model.train(
    ..., accelerator="gpu", devices=-1, strategy="ddp_find_unused_parameters_true"
)
```
- From interactive session (e.g. jupyter notebook):
```python
model.train(
    ...,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_notebook_find_unused_parameters_true",
)
```

3. **Shared memory for data deduplication:**

   By default, shared memory is auto-enabled when DDP is detected (`share_memory=None`). You can explicitly control this behavior:

   ```python
   # Explicitly enable shared memory
   model.train(
       ...,
       accelerator="gpu",
       devices=-1,
       strategy="ddp_find_unused_parameters_true",
       datasplitter_kwargs={"share_memory": True},
   )

   # Explicitly disable shared memory
   model.train(
       ...,
       accelerator="gpu",
       devices=-1,
       strategy="ddp_find_unused_parameters_true",
       datasplitter_kwargs={"share_memory": False},
   )
   ```

   :::{note}
   Shared memory deduplication only applies to dense numpy or scipy sparse `adata.X`. If `adata` is backed (h5ad on disk) or uses dask arrays, shared memory is automatically skipped since those formats already avoid full in-memory copies.
   :::

4. There are a few limitations with the current implementation:
   - During an interactive session, like in a jupyter notebook, we can only train 1 model in multi GPU mode, per session.
   It means that we can't train SCANVI model from SCVI model if the SCVI model was trained in the same notebook. Therefore, need to train and save the SCVI model in another session and load it in the other session. This is a torch lightning caveat.
   - It can't run with early stopping right now (and some models, like totalvi, use early stopping by default), so we disable early stopping once running with DDP. the reason is that the validation loop should be running on one device only and not multiGPU.
Those caveats fixes are still under development as to SCVI-Tools v1.3.
