FROM nvidia/cuda:12.5.0-runtime-ubuntu22.04
FROM python:3.11 AS base

RUN pip install --no-cache-dir uv

RUN uv pip install --system --no-cache torch torchvision torchaudio "jax[cuda12]"

CMD ["/bin/bash"]

FROM base AS build

ENV SCVI_PATH="/usr/local/lib/python3.11/site-packages/scvi-tools"

COPY . ${SCVI_PATH}

ARG DEPENDENCIES=""
RUN uv pip install --system "scvi-tools[${DEPENDENCIES}] @ ${SCVI_PATH}"
