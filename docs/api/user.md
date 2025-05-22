# User

Import scvi-tools as:

```
import scvi
```

```{eval-rst}
.. currentmodule:: scvi

```

## Model

```{eval-rst}
.. currentmodule:: scvi

```

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   model.AUTOZI
   model.CondSCVI
   model.DestVI
   model.LinearSCVI
   model.PEAKVI
   model.SCANVI
   model.SCVI
   model.TOTALVI
   model.MULTIVI
   model.AmortizedLDA
   model.JaxSCVI
```

```{eval-rst}
.. currentmodule:: scvi
```

## External models

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   external.CellAssign
   external.GIMVI
   external.RNAStereoscope
   external.SpatialStereoscope
   external.SOLO
   external.SCAR
   external.Tangram
   external.SCBASSET
   external.ContrastiveVI
   external.POISSONVI
   external.VELOVI
   external.MRVI
   external.METHYLVI
   external.METHYLANVI
   external.Decipher
   external.TOTALANVI
   external.RESOLVI
   external.SysVI
```

## Data loading

`scvi-tools` relies entirely on the [AnnData] format. For convenience, we have included data loaders from the [AnnData] API. [Scanpy] also has [utilities] to load data that are outputted by 10x's Cell Ranger software.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.read_h5ad
   data.read_csv
   data.read_loom
   data.read_text
   data.read_10x_atac
   data.read_10x_multiome

```

## Basic preprocessing

For general single-cell preprocessing, we defer to our friends at [Scanpy], and specifically their preprocessing module ({mod}`scanpy.pp`).

All `scvi-tools` models require raw UMI count data. The count data can be safely stored in an AnnData layer as one of the first steps of a Scanpy single-cell workflow:

```
adata.layers["counts"] = adata.X.copy()
```

Here we maintain a few package specific utilities for feature selection, etc.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.poisson_gene_selection
   data.organize_cite_seq_10x
   data.organize_multiome_anndatas
   data.add_dna_sequence
   data.reads_to_fragments
```

```{eval-rst}
.. currentmodule:: scvi
```

## Model hyperparameter tuning

`scvi-tools` supports automatic model hyperparameter tuning using [Ray Tune].

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   autotune.run_autotune
   autotune.AutotuneExperiment
```

## Model hub

We have a hub for pre-trained `scvi-tools` models that is hosted on [huggingface](https://huggingface.co/models).
Using the functionality that `scvi-tools` provides, users can download pre-trained `scvi-tools` models (and datasets)
from this platform, and model generators can upload their own pre-trained `scvi-tools` models to this platform.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   hub.HubMetadata
   hub.HubModelCardHelper
   hub.HubModel
```

## Model criticism

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   criticism.PosteriorPredictiveCheck
```

## Utilities

Here we maintain miscellaneous general methods.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   model.utils.get_minified_adata_scrna
```

## Configuration

An instance of the {class}`~scvi._settings.ScviConfig` is available as `scvi.settings` and allows configuring scvi-tools.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   _settings.ScviConfig
```

[anndata]: https://anndata.readthedocs.io/en/stable/
[scanpy]: https://scanpy.readthedocs.io/en/stable/index.html
[utilities]: https://scanpy.readthedocs.io/en/stable/api/index.html#reading
[ray tune]: https://docs.ray.io/en/latest/tune/index.html
