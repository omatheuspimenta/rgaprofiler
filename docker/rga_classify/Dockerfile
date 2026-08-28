# rga_classify -- vendored from omatheuspimenta/SugarcaneTranscriptomics
# (code/rgas/, https://github.com/omatheuspimenta/SugarcaneTranscriptomics),
# the reference RGA-classification pipeline this whole tool exists to feed.
#
# This is a straight vendored copy (docker/rga_classify/src/, code + config
# only, no test suite) run exactly as documented upstream -- `uv sync` then
# `uv run python code/rgas/rgas_prediction.py`, no patches, no rewritten
# classification rules. Nothing here is license-gated (MIT-style project,
# pure Python + pandas/numpy/PyYAML/rich, no model weights), so unlike
# deeploc2/signalp6/deeptmhmm this image needs nothing mounted from
# softwares/ at runtime.
#
# Base image bundles both uv and the pinned Python (pyproject.toml requires
# >=3.14) so no separate python3.14 install step is needed -- see
# https://github.com/astral-sh/uv/pkgs/container/uv.
FROM ghcr.io/astral-sh/uv:python3.14-bookworm-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_PREFERENCE=only-managed \
    UV_PYTHON_INSTALL_DIR=/opt/rga_classify/.uv-python \
    UV_CACHE_DIR=/opt/rga_classify/.uv-cache

# procps (for `ps`) is required by Nextflow itself to collect task resource
# metrics inside the container.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/rga_classify
COPY src/pyproject.toml src/uv.lock ./
RUN uv sync --locked --no-install-project

COPY src/code ./code

RUN groupadd --gid 1000 rga_classify \
    && useradd --uid 1000 --gid rga_classify --create-home --shell /bin/bash rga_classify \
    && chmod -R a+rwX /opt/rga_classify

# Nextflow's docker profile runs containers as the host UID, not this image's
# own uid-1000 user -- /opt/rga_classify/.venv (created by `uv sync` above)
# and the outdir the module writes to must stay writable by an arbitrary UID.
USER rga_classify
WORKDIR /data
