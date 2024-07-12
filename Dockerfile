FROM nvidia/cuda:12.5.0-runtime-ubuntu22.04
FROM python:3.11 AS base

RUN pip install --no-cache-dir uv

RUN uv pip install --system --no-cache torch torchvision torchaudio "jax[cuda12]"

CMD ["/bin/bash"]

FROM base AS latest
ARG DEPENDENCIES=""
RUN uv pip install --system --no-cache \
    "scvi-tools[${DEPENDENCIES}] @ git+https://github.com/scverse/scvi-tools"

FROM base AS stable
ARG DEPENDENCIES=""
RUN uv pip install --system --no-cache \
    "scvi-tools[${DEPENDENCIES}]"

FROM base AS branch
ARG SCVI_TOOLS_VERSION=main
ARG DEPENDENCIES=""
RUN uv pip install --system --no-cache \
    "scvi-tools[${DEPENDENCIES}] @ git+https://github.com/scverse/scvi-tools@${SCVI_TOOLS_VERSION}"

FROM base AS semver
ARG SCVI_TOOLS_VERSION
ARG DEPENDENCIES=""
RUN uv pip install --system --no-cache "scvi-tools[${DEPENDENCIES}]==${SCVI_TOOLS_VERSION}"
