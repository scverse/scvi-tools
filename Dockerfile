FROM python:3.13 AS base

ENV LD_LIBRARY_PATH=/usr/local/lib/python3.13/site-packages/nvidia/cu13/lib:/usr/local/nvidia/lib64:/usr/local/nvidia/lib

RUN pip install --no-cache-dir uv

RUN uv pip install --system --no-cache torch torchvision torchaudio

CMD ["/bin/bash"]

FROM base AS build

ENV SCVI_PATH="/usr/local/lib/python3.13/site-packages/scvi-tools"

COPY . ${SCVI_PATH}

ARG DEPENDENCIES=""
RUN uv pip install --system "scvi-tools[${DEPENDENCIES}] @ ${SCVI_PATH}"
