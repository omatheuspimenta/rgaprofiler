# DeepCoil (labstructbioinf/DeepCoil) container
#
# DeepCoil's PyPI release pins a 2020-era ML stack.
# The repo requires python>=3.7,<3.9 -> we use 3.7 to fetch pre-compiled wheels
# for older packages like spacy 2.1.9 and avoid legacy build-from-source errors.
#
# Multi-stage: the builder stage needs a C/C++ toolchain to compile jsonnet
# (an allennlp dependency); the runtime stage only needs the resulting venv
# plus libopenblas's shared library, so the toolchain never ships in the
# final image.

FROM python:3.7-slim-bullseye AS builder

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000

# Build tools are needed to compile jsonnet (allennlp dependency) from source.
# libopenblas speeds up numpy; unzip is needed for the seqvec weights below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        cmake \
        curl \
        unzip \
        libopenblas-dev \
        libssl-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /opt/deepcoil

# numpy + cython pinned first.
RUN pip install --upgrade pip==23.3.2 setuptools==65.5.1 wheel && \
    pip install \
        numpy==1.18.5 \
        "cython<3.0"

# TensorFlow. The original 2.3.0 pin (deepcoil's 2020-era stack) needs CUDA 10.1, which
# doesn't support Ampere-and-later GPUs (compute capability >= 8.0) at all -- not just
# "libraries missing", the CUDA version itself can't target this hardware. 2.11.0 is the
# newest release that still ships cp37 wheels (matching this image's Python 3.7 pin) and
# wants CUDA 11.2/cuDNN 8.1, which does support Ampere. The nvidia-*-cu11 pip packages
# below are TF's own redistributable CUDA/cuDNN runtime libraries -- no system CUDA
# toolkit or alternate base image needed, only LD_LIBRARY_PATH (set in the runtime stage).
RUN pip install tensorflow==2.11.0 \
        nvidia-cudnn-cu11==8.6.0.163 \
        nvidia-cublas-cu11 \
        nvidia-cuda-nvrtc-cu11 \
        nvidia-cuda-runtime-cu11 \
        nvidia-cufft-cu11 \
        nvidia-curand-cu11 \
        nvidia-cusolver-cu11 \
        nvidia-cusparse-cu11

# torch, pinned within deepcoil's torch<2.0,>=1.2 range. deepcoil's internal allennlp/ELMo
# encoder (deepcoil.utils.seqvec.SeqVec) runs on torch, not just TF -- 1.7.1's default
# PyPI wheel predates Ampere GPU support ("CUDA error: no kernel image is available for
# execution on the device" at runtime, confirmed empirically), so bumped to the newest
# version still inside deepcoil's declared upper bound; 1.13.1's default wheel has
# Ampere/Ada kernels built in.
RUN pip install torch==1.13.1

# jsonnet pinned to an older version (2020 era) so allennlp
# doesn't fetch a modern version that requires Python 3.8+ build tools.
RUN pip install jsonnet==0.15.0

# Pin spaCy's C-extensions to their 2019 versions. This prevents
# pip from fetching modern source distributions that require Cython >= 3.1.
RUN pip install \
        "murmurhash==1.0.2" \
        "cymem==2.0.2" \
        "preshed==2.0.1" \
        "blis==0.2.2" \
        "thinc==7.0.8"

# allennlp 0.9.0 and its overrides dependency.
# Modern 'overrides' (v7+) does strict runtime type-checking that breaks
# allennlp 0.9.0. Pinning to a 3.x version fixes the TypeError.
RUN pip install overrides==3.1.0 allennlp==0.9.0

# the remaining direct deepcoil dependencies.
RUN pip install \
        "biopython>=1.77" \
        "pandas==1.1.5" \
        "seaborn==0.11.2" \
        "matplotlib==3.3.4"

# deepcoil itself.
RUN pip install deepcoil==2.0.1 --no-deps

# DeepCoil's ESM/seqvec weights bundle, fetched into the venv's site-packages.
RUN mkdir -p /opt/venv/lib/python3.7/site-packages/deepcoil/weights/seqvec && \
    curl -k -L -o /tmp/seqvec.zip https://lbs.cent.uw.edu.pl/static/files/seqvec.zip && \
    unzip /tmp/seqvec.zip -d /opt/venv/lib/python3.7/site-packages/deepcoil/weights/seqvec && \
    rm /tmp/seqvec.zip

# --- Runtime stage: only what's needed to execute `deepcoil` ---
FROM python:3.7-slim-bullseye AS runtime

# allennlp's dependency chain (via tensorboardX) needs protobuf>=3.20 (its generated
# _pb2.py uses the google.protobuf.internal.builder module added in 3.20), but
# tensorflow==2.11.0's own pre-generated _pb2.py files were compiled pre-3.20 and error
# under protobuf's newer default C++ backend ("Descriptors cannot not be created
# directly"). No single protobuf version satisfies both. Forcing the pure-Python
# implementation (protobuf's own documented workaround for exactly this situation)
# resolves it without pinning protobuf down and breaking tensorboardX instead.
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# libopenblas-dev's runtime .so is needed by numpy/tensorflow at execution time.
# procps (for `ps`) is required by Nextflow itself to collect task resource metrics
# inside the container — without it, task execution fails with "Command 'ps'
# required by nextflow to collect task metrics cannot be found".
RUN apt-get update && apt-get install -y --no-install-recommends \
        libopenblas-dev \
        procps \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 1000 deepcoil \
    && useradd --uid 1000 --gid deepcoil --create-home --shell /bin/bash deepcoil

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

# TensorFlow's GPU device discovery needs these pip-installed CUDA/cuDNN runtime
# libraries (see the tensorflow install step above) on the loader path.
ENV LD_LIBRARY_PATH="/opt/venv/lib/python3.7/site-packages/nvidia/cudnn/lib:/opt/venv/lib/python3.7/site-packages/nvidia/cublas/lib:/opt/venv/lib/python3.7/site-packages/nvidia/cuda_runtime/lib:/opt/venv/lib/python3.7/site-packages/nvidia/cufft/lib:/opt/venv/lib/python3.7/site-packages/nvidia/curand/lib:/opt/venv/lib/python3.7/site-packages/nvidia/cusolver/lib:/opt/venv/lib/python3.7/site-packages/nvidia/cusparse/lib"

USER deepcoil
WORKDIR /data
