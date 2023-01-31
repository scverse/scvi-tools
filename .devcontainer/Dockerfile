# syntax=docker/dockerfile:1
FROM python:3.10
WORKDIR /scvi-tools
COPY . .

RUN pip install --upgrade pip
RUN pip install "https://storage.googleapis.com/jax-releases/nocuda/jaxlib-0.4.1-cp310-cp310-manylinux2014_x86_64.whl"
RUN pip install -e ".[dev,docs]"
