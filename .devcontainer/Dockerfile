# syntax=docker/dockerfile:1
FROM python:3.10
WORKDIR /scvi-tools
COPY . .

RUN pip install --upgrade pip
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
RUN bash Mambaforge-$(uname)-$(uname -m).sh
RUN pip install -e ".[dev,docs,tutorials]"
