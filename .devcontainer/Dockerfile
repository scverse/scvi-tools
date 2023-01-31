# syntax=docker/dockerfile:1
FROM python:3.10
WORKDIR /scvi-tools

COPY .. .

RUN pip install --upgrade pip
RUN pip install -e ".[dev,docs]"
