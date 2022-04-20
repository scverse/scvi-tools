# Installation

## Prerequisites

scvi-tools can be installed via conda or pip. If you don't know which to choose, we recommend conda for beginner users.

### conda prerequisites

1. Install Conda. We typically use the [Miniconda] Python distribution. Use Python version >=3.7.

2. Create a new conda environment:

   ```
   conda create -n scvi-env python=3.9
   ```

3. Activate your environment:

   ```
   conda activate scvi-env
   ```

### pip prerequisites

1. Install [Python], we prefer the [pyenv](https://github.com/pyenv/pyenv/) version management system, along with [pyenv-virtualenv](https://github.com/pyenv/pyenv-virtualenv/).
2. Install [PyTorch]. If you have an Nvidia GPU, be sure to install a version of PyTorch that supports it -- scvi-tools runs much faster with a discrete GPU.

:::{note}
Installing scvi-tools on a Mac with Apple Silicon is only possible using a native version of Python. A native version of Python can be installed with an Apple Silicon version of miniconda (which can be installed from a native version of homebrew). This is due to an scvi-tools dependency on jax, which cannot be run via Rosetta. 
:::

## Conda

```
conda install scvi-tools -c conda-forge
```

## Pip

```
pip install scvi-tools
```

Through pip with packages to run notebooks. This installs scanpy, etc.:

```
pip install scvi-tools[tutorials]
```

Nightly version - clone this repo and run:

```
pip install .
```

## Development

For development - clone this repo and run:

```
pip install -e ".[dev,docs]"
```

## R

scvi-tools can be called from R via Reticulate.

This is only recommended for basic functionality (getting the latent space, normalized expression, differential expression). For more involved analyses with scvi-tools, we highly recommend using it from Python.

The easiest way to install scvi-tools for R is via conda.

1. Install Conda Prerequisites (see above).

2. Install Reticulate:

   ```
   install.packages("reticulate")
   ```

3. Then in your R code:

   ```
   library(reticulate)
   use_condaenv("scvi-env", required=TRUE)
   ```

[miniconda]: https://conda.io/miniconda.html
[python]: https://www.python.org/downloads/
[pytorch]: http://pytorch.org
[reticulate]: https://rstudio.github.io/reticulate/
