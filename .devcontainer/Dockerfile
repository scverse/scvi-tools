# syntax=docker/dockerfile:1
FROM condaforge/mambaforge:22.9.0-3
WORKDIR /scvi-tools
COPY . .

RUN mamba init
RUN mamba install python=${PYTHON_VERSION} --yes
RUN mamba install scvi-tools --yes
RUN pip install -e ".[dev,docs,autotune,hub]"
