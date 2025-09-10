# Installation

## Quick install

scvi-tools can be installed via `conda` or `pip`.
We recommend installing into a fresh virtual environment to avoid conflicts with other packages
and compatibility issues.

For the basic CPU version run:

```bash
pip install -U scvi-tools
```
or
```bash
conda install scvi-tools -c conda-forge
```

In order to install scvi-tools with Nvidia GPU CUDA support, for Linux Systems
(such as Ubuntu or RedHat) use:

```bash
pip install -U scvi-tools[cuda]
```
And for Apple Silicon metal (MPS) support:
```bash
pip install -U scvi-tools[metal]
```

Don't know how to get started with virtual environments or `conda`/`pip`? Check out the
[prerequisites](#prerequisites) section.

## Prerequisites

### Virtual environment

A virtual environment can be created with either `conda` or `venv`. We recommend using a fresh `conda` environment.
We currently support Python 3.11 - 3.13.

For `conda`, we recommend using the [Miniforge](https://github.com/conda-forge/miniforge) or
[Mamba](https://mamba.readthedocs.io/en/latest/) distribution, which are generally lighter & faster
than the official distribution and comes with conda-forge as the default channel
(where scvi-tools is hosted).

```bash
conda create -n scvi-env python=3.13  # any python 3.11 to 3.13
conda activate scvi-env
```

For `venv`, we recommend using [uv](https://github.com/astral-sh/uv), which is a high-performance
Python package manager and installer written in Rust.

```bash
pip install -U uv
uv venv .scvi-env
source .scvi-env/bin/activate  # for macOS and Linux
.scvi-env\Scripts\activate  # for Windows
```

### GPU support with PyTorch and JAX

scvi-tools depends on PyTorch for accelerated computing (and optionally on Jax). If you don't plan
on using an accelerated device, we recommend installing scvi-tools directly and letting these
dependencies be installed automatically by your package manager of choice.

If you plan on taking advantage of an accelerated device (e.g. Nvidia GPU or Apple Silicon),
which is likely, scvi-tools supports it and you should install with the GPU support dependency of scvi-tools.

However, there might be cases where the GPU HW is not supporting the latest installation of PyTorch and Jax.
In this case we recommend installing PyTorch and JAX _before_ installing scvi-tools.
Please follow the respective installation instructions for [PyTorch](https://pytorch.org/get-started/locally/) and
[JAX](https://jax.readthedocs.io/en/latest/installation.html) compatible with your system and device type.

## Optional dependencies

scvi-tools is installed in its lightest form by default.
It has many optional dependencies, which expand its capabilities:

- _autotune_ - in order to run scvi.autotune
- _hub_ - in order to use scvi.hub
- _regseq_ - in order to run scvi.data.add_dna_sequence
- _file_sharing_ - for convenient files sharing
- _parallel_ - for parallelization engine
- _interpretability_ - for supervised models interpretability
- _dataloaders_ - for custom dataloaders use
- _jax_ - for Jax support
- _tests_ - in order to be able to perform tests
- _editing_ - for code editing
- _dev_ - for development purposes
- _cuda_ - for linux based OS GPU support
- _metal_ - for Apple Silicon metal (MPS) support
- _docsbuild_ - in order to create docs

The easiest way to install these is with `pip`.
In order to install capability X run: _pip install scvi-tools[X]_

You can install several capabilities together, e.g:
To install scvi-tools with JAX support for GPU on Ubuntu: _pip install scvi-tools[cuda,jax]_

To install all tutorial dependencies:

```bash
pip install -U scvi-tools[tutorials]
```

To install all optional dependencies (_e.g._ jax support, custom dataloaders, autotune, criticism, model hub):


```bash
pip install -U scvi-tools[optional]
```

To install development dependencies, including `pre-commit` and testing dependencies:

```bash
pip install -U scvi-tools[dev]
```

And to install all capabilites of scvi-tools, simply run:

```bash
pip install -U scvi-tools[all]
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
