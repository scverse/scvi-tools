# Train SCVI model with custom dataloaders

:::{note}
This page is under construction.
:::

In SCVI, custom dataloaders allow you to create a tailored data pipeline that can handle unique formats or complex datasets not covered by the default loaders. A custom dataloader can be useful when you have a specific structure for your data or need to preprocess it in a particular way before feeding it into the model, in order to gain some advantage.

For example, we offer custom dataloaders that do not necessarily store the data on memory while training, thus enable us to expand the sizes of dataset that we can train our models based on while not being limited by the amount of memory that we have.
Without dataloader large data can be on disk but inefficient. We increase efficiency for data on disk and enable data on cloud storage.

In SCVI, we work with several collaborators in order to construct efficient custom dataloaders:
1. [lamin.ai]("https://lamin.ai/") custom dataloader is based on MappedCollectionDataModule and can run a collection of adata based on lamindb backend.

LamindDB is a key-value store designed specifically for machine learning, particularly focusing on making it easier to manage large-scale training datasets. It is optimized for storing and querying tabular data, providing fast read and write access. LamindDB’s design allows for efficient handling of large datasets with a focus on machine learning workflows, such as those used in SCVI.

Pros:

- Fast Retrieval: LamindDB provides fast access to tabular data, which can be beneficial for training large machine learning models like SCVI on single-cell RNA-seq data.
- Simplicity: Since LamindDB is designed with machine learning in mind, it may offer a simpler interface for managing and querying tabular data compared to TileDB.
- Optimized for ML Workflows: If your dataset is structured as tables (rows and columns), LamindDB’s format aligns well with SCVI's expectations, potentially reducing the need for complex transformations.

```python
os.system("lamin init --storage ./test-registries")
import lamindb as ln

ln.setup.init(name="lamindb_instance_name", storage=save_path)

# a test for mapped collection
collection = ln.Collection.get(name="covid_normal_lung")
artifacts = collection.artifacts.all()
artifacts.df()

datamodule = MappedCollectionDataModule(
    collection, batch_key="assay", batch_size=1024, join="inner"
)
model = scvi.model.SCVI(adata=None, registry=datamodule.registry)
...
```
LamindDB may not be as efficient or flexible as TileDB for handling complex multi-dimensional data

2. [CZI]("https://chanzuckerberg.com/") based [tiledb]("https://tiledb.com/") custom dataloader is based on CensusSCVIDataModule and can run a large multi-dimensional datasets that are stored in TileDB’s format.

TileDB is a general-purpose, multi-dimensional array storage engine designed for high-performance, scalable data access. It supports various data types, including dense and sparse arrays, and is optimized for handling large datasets efficiently. TileDB’s strength lies in its ability to store and query data across multiple dimensions and scale efficiently with large volumes of data.

Pros:

- Efficient Data Access: With TileDB, you can query and access specific subsets of your data without needing to load the entire dataset into memory, improving performance.
Scalability: Handles large datasets that exceed your system's memory capacity, making it ideal for single-cell RNA-seq datasets with millions of cells.
- Flexibility: Supports multi-dimensional data storage, which can be useful for storing gene expression data across multiple conditions or time points.

```python
import cellxgene_census
import tiledbsoma as soma
from cellxgene_census.experimental.ml import experiment_dataloader
from cellxgene_census.experimental.ml.datamodule import CensusSCVIDataModule
import numpy as np

# this test checks the local custom dataloder made by CZI and run several tests with it
census = cellxgene_census.open_soma(census_version="stable")

experiment_name = "mus_musculus"
obs_value_filter = (
    'is_primary_data == True and tissue_general in ["kidney"] and nnz >= 3000'
)

hv_idx = np.arange(100)  # just ot make it smaller and faster for debug

# this is CZI part to be taken once all is ready
batch_keys = ["dataset_id", "assay", "suspension_type", "donor_id"]
datamodule = CensusSCVIDataModule(
    census["census_data"][experiment_name],
    measurement_name="RNA",
    X_name="raw",
    obs_query=soma.AxisQuery(value_filter=obs_value_filter),
    var_query=soma.AxisQuery(coords=(list(hv_idx),)),
    batch_size=1024,
    shuffle=True,
    batch_keys=batch_keys,
    dataloader_kwargs={"num_workers": 0, "persistent_workers": False},
)


# basicaly we should mimiC everything below to any model census in scvi
adata_orig = synthetic_iid()
scvi.model.SCVI.setup_anndata(adata_orig, batch_key="batch")
model = scvi.model.SCVI(adata_orig)
...
```
Key Differences between them in terms of Custom Dataloaders:
1. Data Format:

- TileDB is more flexible for handling complex, multi-dimensional datasets (like sparse matrices with many dimensions), making it ideal for large-scale genomics data that don’t fit neatly into a simple table format.
- LamindDB is optimized for tabular datasets (rows and columns), so it's better suited for single-cell RNA-seq data that is already structured in this way.

2. Efficiency:

- TileDB is highly optimized for large, sparse datasets and supports sophisticated indexing, making it ideal for querying and accessing large portions of a dataset efficiently, even across multiple dimensions.
- LamindDB is optimized for fast tabular data access and may be easier to set up and use in machine learning workflows when dealing with simple table-based data.

3. Complexity:

- TileDB requires more effort to integrate into a custom dataloader, as it’s more complex and can involve handling multi-dimensional data storage.
- LamindDB offers a simpler integration for custom dataloaders when dealing with tabular data and might be more intuitive for machine learning tasks.

When to Use Each:
- TileDB is ideal when your data is large, multi-dimensional, and potentially sparse, requiring advanced indexing and querying capabilities.
- LamindDB is a better choice when you're working with large, structured tabular data and need fast, efficient access for machine learning models like SCVI.

Writing custom dataloaders requires a good understanding of PyTorch’s DataLoader class and how to integrate it with SCVI, which may be difficult for beginners.
It will also requite maintenance: If the data format or preprocessing needs change, you’ll have to modify and maintain the custom dataloader code, But it can be a greate addition to the model pipeline, in terms of runtime and how much data we can digest.

:::{note}
As for SCVI-Tools v1.3.0 Custom Dataloaders are experimental.
:::
