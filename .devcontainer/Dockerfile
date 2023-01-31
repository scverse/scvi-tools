# syntax=docker/dockerfile:1
FROM condaforge/mambaforge:latest
WORKDIR /scvi-tools
COPY . .

RUN mamba init
RUN mamba install python=3.10
RUN mamba install scvi-tools
RUN pip install -e ".[dev,docs,tutorials]"
