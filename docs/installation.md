# Installation

## Quick install

scvi-tools can be installed via `conda` or `pip`. We recommend installing into a virtual
environment to avoid conflicts with other packages.

```bash
conda install scvi-tools -c conda-forge
```

or

```bash
pip install -U scvi-tools
```

Don't know how to get started with virtual environments or `conda`/`pip`? Check out the
[prerequisites](#prerequisites) section.

## Prerequisites

### Virtual environment

A virtual environment can be created with either `conda` or `venv`. We recommend using `conda`. We
currently support Python 3.10 - 3.13.

For `conda`, we recommend using the [Miniforge](https://github.com/conda-forge/miniforge)
distribution, which is generally faster than the official distribution and comes with conda-forge
as the default channel (where scvi-tools is hosted).

```bash
conda create -n scvi-env python=3.12  # any python 3.10 to 3.13
conda activate scvi-env
```

For `venv`, we recommend using [uv](https://github.com/astral-sh/uv).

```bash
pip install -U uv
uv venv .scvi-env
source .scvi-env/bin/activate  # for macOS and Linux
.scvi-env\Scripts\activate  # for Windows
```

### PyTorch and JAX

scvi-tools depends on PyTorch and JAX for accelerated computing. If you don't plan on using
an accelerated device, we recommend installing scvi-tools directly and letting these dependencies
be installed automatically by your package manager of choice.

If you plan on taking advantage of an accelerated device (e.g. Nvidia GPU or Apple Silicon), scvi-tools supports it.
In order to install scvi-tools with Nvidia GPU CUDA support use:
```bash
pip install -U scvi-tools[cuda]
```
And for Apple Silicon metal (MPS) support:
```bash
pip install -U scvi-tools[metal]
```

However, there might be cases where the GPU HW is not supporting the latest installation of PyTorch and Jax.
In this case we recommend installing PyTorch and JAX _before_ installing scvi-tools.
Please follow the respective installation instructions for [PyTorch](https://pytorch.org/get-started/locally/) and
[JAX](https://jax.readthedocs.io/en/latest/installation.html) compatible with your system and device type.

## Optional dependencies

scvi-tools has optional dependencies. The easiest way to install these is with `pip`.

To install all tutorial dependencies:

```bash
pip install -U scvi-tools[tutorials]
```

To install all optional dependencies (_e.g._ autotune, criticism, model hub):

```bash
pip install -U scvi-tools[optional]
```

To install development dependencies, including `pre-commit` and testing dependencies:

```bash
pip install -U scvi-tools[dev]
```

## Docker

If you plan on running scvi-tools in a containerized environment, we provide various Docker
[images](https://github.com/scverse/scvi-tools/pkgs/container/scvi-tools) hosted on GHCR.

## R

scvi-tools can be called from R via Reticulate.

This is only recommended for basic functionality (getting the latent space, normalized expression,
differential expression). For more involved analyses with scvi-tools, we highly recommend using it
from Python.

The easiest way to install scvi-tools for R is via conda.

1. Install conda prerequisites.

2. Install R and reticulate in the conda environment:

    ```bash
    conda install -c conda-forge r-base r-essentials r-reticulate
    ```

3. Then in your R code:

    ```R
    library(reticulate)
    ```
