# Train SCVI model on multiple GPUs


scvi-tools v1.3.0 supports training on a **multi GPU system**, which can significantly speed up training and allow you to handle larger datasets. It is supported only on NVidia GPUs and using Distributed Data Parallel (DDP) within CUDA.

CUDA (Compute Unified Device Architecture) is used to enable high-performance parallel computing on NVIDIA GPUs. DDP spreads the data across multiple GPUs to improve training speed or increase batch size.

Using **CUDA** with **DDP** for multi-GPU training can significantly boost the performance and scalability of your model training, especially with large datasets and models. However, it also comes with its own set of challenges, such as the need for proper configuration, higher hardware costs, and potential software compatibility issues.

## Some of the benfits of using MultiGPU training

### 1. **Faster Training**
   - **Parallel Processing:** Using multiple GPUs allows you to parallelize the training process, meaning different GPUs can work on different batches of data simultaneously, leading to faster training times. However, parallelization requires some overhead and training on 2 GPUs might not provide significant acceleration.
   - **Scalability:** As you add more GPUs to the system, the training speed improves, further reducing the time it takes to train large models.
   - **Larger Models:** With multiple GPUs, you can train models that are too large to fit into a single GPU’s memory. This component is key for recent foundation models, while models within scvi-tools are small enough to not benefit from splitting a model across different GPU's and is not supported currently.

See the following figure which shows the training time of a PBMC data, SCVI model with and without the use of multiGPU vs different cell numbers in data:
We can see the advantage with larger data, while for the small data, theres no much advantage due to multiGPU overhead needed.

:::{figure} /\_static/multigpu.png
:align: center
:alt: MultiGPU train time
:class: img-fluid
:::

### 2. **Efficiency with DDP**
   - **Data Parallelism:** DDP splits your data into smaller chunks, distributes them across multiple GPUs, and then aggregates the results. This allows for better use of available hardware resources.
   - **Synchronization:** DDP efficiently synchronizes gradients and updates between GPUs after each batch, ensuring that each GPU works in sync with the others without much manual intervention.
   - **Automatic Handling of Distributed Work:** DDP abstracts much of the complexity of managing multiple GPUs, so you don’t have to manually manage how data is split or how results are aggregated.


## Using MultiGPU training in SCVI-Tools

1. Currently, multiGPU training is supported for the following models:
   - {class}`~scvi.model.SCVI`
   - {class}`~scvi.model.SCANVI`,
   - {class}`~scvi.model.CondSCVI`
   - {class}`~scvi.model.LinearSCVI`
   - {class}`~scvi.model.TOTALVI`,
   - {class}`~scvi.model.MULTIVI`
   - {class}`~scvi.model.PEAKVI`


2. How to use it: during the model train command, we need to use the strategy parameter and its value depends on whether we are running from command line/script or from an interactive session like jupyter notebook/colab:
- From non-interactive session:
```python
model.train(
    ...,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true"
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

3. There are a few caveats with the current implementation:
   - During interactive session, like in jupyter notebook, we can only train one model in multi GPU mode, per session.
   It means that we can't train SCANVI model from SCVI model, if the SCVI model was trained in the same notebook. Therefore, for example we need to train and save the SCVI model and load it in a new session for training SCANVI on multiple GPUs. This is a pytorch Lightning caveat. Contributions are welcome to fix this behavior.
   - It doesn't support early stopping currently (e.g. totalVI uses early stopping by default). In order to support training, we disable early stopping once running with DDP. The validation loop needs to be executed on a single device only and not using multiple GPUs.
