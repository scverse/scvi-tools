# syntax=docker/dockerfile:1
FROM condaforge/mambaforge:latest
WORKDIR /scvi-tools
COPY . .

RUN mamba create -n scvi-tools-dev python=3.10
RUN mamba init
RUN mamba activate scvi-tools-dev
RUN mamba install scvi-tools
RUN pip install -e ".[dev,docs,tutorials]"
