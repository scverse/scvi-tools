# Installation

## Prerequisites

scvi-tools can be installed via conda or pip. If you don't know which to choose, we recommend conda.

### conda

1. Install conda. We typically use the [mambaforge] distribution. Use python>=3.9, and consider using mamba instead of conda. Mamba is a drop-in replacement for conda that is significantly more efficient.

2. Create a new conda environment:

    ```
    conda create -n scvi-env python=3.9
    ```

3. Activate your environment:

    ```
    conda activate scvi-env
    ```

### pip

1. If using conda/mamba, then just run `conda install -c anaconda pip` and skip this section.
2. Install [Python], we prefer the [pyenv](https://github.com/pyenv/pyenv/) version management system, along with [pyenv-virtualenv](https://github.com/pyenv/pyenv-virtualenv/).
3. Install [PyTorch] and [jax]. If you have an Nvidia GPU, be sure to install versions of PyTorch and jax that support it -- scvi-tools runs much faster with a discrete GPU.

### Apple silicon

Installing scvi-tools on a Mac with Apple Silicon is only possible using a native version of python. A native version of python can be installed with an Apple Silicon version of mambaforge (which can be installed from a native version of homebrew via `brew install --cask mambaforge`). This is due to an scvi-tools dependency on jax, which cannot be run via Rosetta.

### Windows

After setting up a virtual environment with conda/mamba/pyenv, please install the [community-built version of jax](https://github.com/cloudhan/jax-windows-builder) before installing scvi-tools.

```
pip install "jax[cpu]" -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver
```

### GPU

All scvi-tools models will be faster when accelerated with a GPU. Before installing scvi-tools, you can install GPU versions of PyTorch and jax using conda as follows:

```
conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
conda install jax jaxlib -c conda-forge
```

Please go to the respective package website for more information on how to install with pip.

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
pip install "scvi-tools[tutorials]"
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

1. Install conda prerequisites.

2. Install R and reticulate in the conda environment:

    ```
    conda install -c conda-forge r-base r-essentials r-reticulate
    ```

3. Then in your R code:

    ```
    library(reticulate)
    ```

[mambaforge]: https://github.com/conda-forge/miniforge
[python]: https://www.python.org/downloads/
[pytorch]: http://pytorch.org
[jax]: https://jax.readthedocs.io/en/latest/
[reticulate]: https://rstudio.github.io/reticulate/
